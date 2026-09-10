import argparse
import asyncio
import json
import re
from os import environ
from pathlib import Path

import aiofiles
import aiohttp
import dotenv
import pandas as pd
import requests
from loguru import logger

from src.config import get_settings

_MAX_TOKEN_RETRIES = 3
_HTTP_OK = 200
_HTTP_UNAUTHORIZED = 401
_REQUEST_TIMEOUT = 30


def get_token() -> str | None:
    """Get authentication token from GeneStory API using credentials from .env file.

    Returns:
        Authentication token if successful, None otherwise
    """
    # Load environment variables
    dotenv.load_dotenv()

    # Get credentials from environment
    username = environ.get("PORTAL_USERNAME")
    password = environ.get("PORTAL_PASSWORD")

    if not username or not password:
        logger.error("Username or password not found in .env file")
        return None

    url = "https://genestory.ai/gt-authen/users/signin"

    payload = json.dumps(
        {
            "type": "phone_number",
            "username": username,
            "region_code": "+84",
            "password": password,
        },
    )

    headers = {
        "accept": "application/json, text/plain, */*",
        "accept-language": "en-US,en;q=0.9,vi-VN;q=0.8,vi;q=0.7",
        "authorization": "Bearer undefined",
        "content-type": "application/json",
        "origin": "https://genestory.ai",
        "priority": "u=1, i",
        "referer": "https://genestory.ai/cs",
        "sec-ch-ua": '"Google Chrome";v="137", "Chromium";v="137", "Not/A)Brand";v="24"',
        "sec-ch-ua-mobile": "?0",
        "sec-ch-ua-platform": '"Linux"',
        "sec-fetch-dest": "empty",
        "sec-fetch-mode": "cors",
        "sec-fetch-site": "same-origin",
        "user-agent": (
            "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/137.0.0.0 Safari/537.36"
        ),
        "Cookie": (
            "_gid=GA1.2.341526983.1748833596; wp-wpml_current_language=vi; "
            "_ga_VX1389DHCK=GS2.1.s1749175918$o13$g1$t1749175918$j60$l0$h0; "
            "_ga=GA1.2.1901061896.1748246285; _gat_gtag_UA_229650118_3=1"
        ),
    }

    for retry in range(_MAX_TOKEN_RETRIES):
        try:
            logger.info(f"Attempting to get token (attempt {retry + 1}/{_MAX_TOKEN_RETRIES})")
            # Send the POST request
            response = requests.post(url, headers=headers, data=payload, timeout=_REQUEST_TIMEOUT)
        except requests.exceptions.RequestException as e:
            logger.error(f"Request failed (attempt {retry + 1}): {e}")
            continue

        # Check if the response is successful
        if response.status_code != _HTTP_OK:
            logger.error(f"Authentication failed - Status: {response.status_code}, Response: {response.text}")
            return None

        response_data = response.json()
        if "data" in response_data and "token" in response_data["data"]:
            token = response_data["data"]["token"]
            logger.success("Successfully obtained authentication token")
            return f"Bearer {token}"
        logger.error(f"Unexpected response format: {response_data}")
        return None

    logger.error(f"Failed to get token after {_MAX_TOKEN_RETRIES} attempts")
    return None


def _read_sample_dataframe(file_path: str, file_path_obj: Path) -> pd.DataFrame:
    """Read a sample manifest Excel or TSV file into a DataFrame."""
    if file_path_obj.suffix.lower() in [".xlsx", ".xls"]:
        # Read Excel file with string data types to preserve leading zeros
        df = pd.read_excel(file_path, dtype=str)

        # Find the relevant columns (adjust column names based on your Excel structure)
        lid_col = df.filter(like="LID").columns[0]
        barcode_col = df.filter(like="BARCODE").columns[0]

        # Extract and format the data - keep as strings to preserve leading zeros
        return pd.DataFrame({"lid": df[lid_col], "barcode": df[barcode_col]})

    if file_path_obj.suffix.lower() == ".tsv":
        # Read TSV file - assume first column is LID, second is barcode
        df = pd.read_csv(
            file_path,
            sep="\t",
            dtype=str,
            header=None,
            names=["lid", "barcode"],
        )
        return df.copy()

    msg = f"Unsupported file format: {file_path_obj.suffix}. Only .xlsx, .xls, and .tsv are supported."
    raise ValueError(msg)


class ReportGenerator:
    """PDF report generator for mtDNA samples"""

    def __init__(self) -> None:
        """Initialize the ReportGenerator"""
        self.auth_token: str | None = None

    def extract_sample_info(
        self,
        file_path: str,
        output_path: str,
        barcode_len: int | None = None,
    ) -> list[tuple[str, str]]:
        """
        Extract LID and barcode from Excel or TSV file and save to TSV.
        Expected Excel format should have columns containing LID and barcode information.
        Expected TSV format: first column is LID, second is barcode.
        Only keeps non-empty, non-control LIDs, filtering out PC/NTC control samples.

        Args:
            file_path: Path to the Excel or TSV file
            output_path: Path to output TSV file
            barcode_len: Target length for barcodes (if None, will be determined from data)

        Returns:
            List of tuples containing (lid, barcode) pairs
        """
        try:
            file_path_obj = Path(file_path)
            result_df = _read_sample_dataframe(file_path, file_path_obj)

            logger.info(f"Read {len(result_df)} rows from {file_path}")

            # Filter out empty LIDs and control samples (PC/NTC prefix)
            result_df["lid_upper"] = result_df["lid"].fillna("").str.strip().str.upper()
            result_df = result_df[
                result_df["lid_upper"].ne("")
                & ~result_df["lid_upper"].str.startswith("PC")
                & ~result_df["lid_upper"].str.startswith("NTC")
            ].drop("lid_upper", axis=1)

            # Ensure barcodes have leading zeros
            # First check if barcode can be converted to integer - also handling NaN values
            numeric_mask = result_df["barcode"].notna() & result_df["barcode"].str.match(r"^\d+$")  # pyright: ignore[reportAttributeAccessIssue]
            if numeric_mask.any():
                # If barcode_len not specified, determine it from the longest existing barcode
                if barcode_len is None:
                    # Find the max length of numeric barcodes
                    max_len = result_df.loc[numeric_mask, "barcode"].str.len().max()
                    barcode_len = max_len

                # For numeric barcodes, ensure they're padded with leading zeros
                result_df.loc[numeric_mask, "barcode"] = (
                    result_df.loc[numeric_mask, "barcode"].astype(str).str.zfill(barcode_len)
                )

            # Save to TSV - use string formatting to preserve leading zeros
            result_df.to_csv(output_path, sep="\t", index=False, header=False)
            logger.info(f"Successfully extracted sample information to {output_path}")
            logger.info(f"Total samples after filtering: {len(result_df)}")
            logger.info(f"Barcode length used for padding: {barcode_len}")

            # Return the list of lid, barcode tuples
            return [(row.lid, row.barcode) for _, row in result_df.iterrows()]

        except Exception as e:
            logger.error(f"Error processing file: {e}")
            raise

    async def _request_pdf_link(
        self,
        session: aiohttp.ClientSession,
        lid: str,
        barcode: str,
        json_path: Path,
        auth_token: str,
    ) -> str:
        """Send a single PDF-link request; raise on auth or HTTP failure."""
        headers = {"Authorization": f"{auth_token}"}

        # Create form data
        form_data = aiohttp.FormData()
        form_data.add_field("kit_id", barcode)

        # Read and add the JSON file
        async with aiofiles.open(json_path, "rb") as f:
            json_data = await f.read()
            form_data.add_field(
                "loci_file",
                json_data,
                filename=f"{lid}.json",
                content_type="application/json",
            )

        async with session.post(
            get_settings().parameters.REPORT_BASE_URL,
            headers=headers,
            data=form_data,
        ) as response:
            if response.status == _HTTP_UNAUTHORIZED:
                logger.error(f"Authentication failed for {lid}. Please check your auth token.")
                msg = "Authentication failed. Please check your auth token."
                raise aiohttp.ClientResponseError(
                    request_info=response.request_info,
                    history=response.history,
                    status=_HTTP_UNAUTHORIZED,
                    message=msg,
                )
            response.raise_for_status()
            # Response wraps the URL in JSON-array-ish brackets/quotes; strip those.
            pdf_link = await response.text()
            return re.sub(r"^[\[\" \]]+|[\[\" \]]+$", "", pdf_link)

    async def get_pdf_link(
        self,
        session: aiohttp.ClientSession,
        lid: str,
        barcode: str,
        json_path: Path,
        auth_token: str,
    ) -> str:
        """Get PDF link from the API."""
        retries = get_settings().parameters.REPORT_RETRIES
        for attempt in range(retries):
            try:
                pdf_link = await self._request_pdf_link(session, lid, barcode, json_path, auth_token)
            except aiohttp.ClientResponseError as e:
                if e.status == _HTTP_UNAUTHORIZED:  # If unauthorized, no point in retrying
                    raise
                if attempt == retries - 1:
                    raise
                logger.warning(f"Attempt {attempt + 1} failed for {lid}: {e}")
                await asyncio.sleep(get_settings().parameters.REPORT_RETRY_DELAY)
            except (aiohttp.ClientError, OSError) as e:
                if attempt == retries - 1:
                    raise
                logger.warning(f"Attempt {attempt + 1} failed for {lid}: {e}")
                await asyncio.sleep(get_settings().parameters.REPORT_RETRY_DELAY)
            else:
                return pdf_link

        msg = f"Failed to get PDF link for {lid} after {retries} attempts"
        raise RuntimeError(msg)

    async def download_pdf(self, session: aiohttp.ClientSession, pdf_link: str, output_path: Path) -> None:
        """Download PDF file from the given link."""
        try:
            async with session.get(pdf_link) as response:
                response.raise_for_status()

                async with aiofiles.open(output_path, "wb") as f:
                    await f.write(await response.read())

                logger.success(f"PDF download completed: {output_path}")
        except Exception as e:
            logger.error(f"Error downloading PDF from {pdf_link}: {e}")
            raise

    def get_existing_pdfs(self, report_dir: Path) -> set[str]:
        """Get a set of sample IDs that already have PDF files."""
        existing_pdfs = set()
        if report_dir.exists():
            for pdf_file in report_dir.glob("*.pdf"):
                # Get the sample ID from the filename (remove .pdf extension)
                sample_id = pdf_file.stem
                existing_pdfs.add(sample_id)

            if existing_pdfs:
                logger.debug(f"Found {len(existing_pdfs)} existing PDFs in {report_dir}")
        else:
            logger.debug(f"Report directory {report_dir} does not exist yet")

        return existing_pdfs

    async def process_sample(
        self,
        session: aiohttp.ClientSession,
        lid: str,
        barcode: str,
        results_dir: Path,
        report_dir: Path,
        auth_token: str,
        existing_pdfs: set[str],
    ) -> bool:
        """Process a single sample by getting PDF link and downloading the PDF.
        Returns True if a new PDF was generated, False if skipped or failed.
        """
        # Check if PDF already exists for this sample
        if lid in existing_pdfs:
            logger.info(f"PDF already exists for {lid}, skipping")
            return False

        try:
            json_path = results_dir / lid / f"{lid}.json"
            if not json_path.exists():
                logger.error(f"JSON file not found for {lid}: {json_path}")
                return False
            pdf_link = await self.get_pdf_link(session, lid, barcode, json_path, auth_token)
            logger.info(f"Got PDF link for {lid}: {pdf_link}")

            await self.download_pdf(session, pdf_link, report_dir / f"{lid}.pdf")
            logger.info(f"Downloaded PDF for {lid}")
        except Exception as e:  # noqa: BLE001
            logger.error(f"Processing sample {lid}: {e}")
            return False
        else:
            return True

    def find_data_file(self, batch_id: str, metadata_dir: Path) -> Path | None:
        """Find the Excel or TSV file for a given batch ID in the metadata directory."""
        # Try exact match first (replacing .tsv with .xlsx)
        excel_pattern = f"{batch_id.replace('.tsv', '')}.xlsx"
        data_files = list(metadata_dir.glob(excel_pattern))

        # If no Excel file found, try TSV file
        if not data_files:
            tsv_pattern = f"{batch_id.replace('.tsv', '')}.tsv"
            data_files = list(metadata_dir.glob(tsv_pattern))

        if not data_files:
            logger.warning(f"No Excel or TSV files found for batch {batch_id} in {metadata_dir}")
            return None

        logger.debug(f"Found data file: {data_files[0]}")
        return data_files[0]

    def _resolve_data_dir(self, data_dir: str | None, base_dir: Path) -> Path:
        """Resolve the data directory: explicit arg > DATA_DIR env > base_dir/data."""
        if data_dir:
            return Path(data_dir)
        env_data_dir = environ.get("DATA_DIR")
        if env_data_dir:
            return Path(env_data_dir)
        return base_dir / "data"

    def _resolve_results_dir(self, results_dir: str | None, base_dir: Path) -> Path:
        """Resolve the results directory: explicit arg > RESULTS_DIR env > base_dir/results."""
        if results_dir:
            return Path(results_dir)
        env_results_dir = environ.get("RESULTS_DIR")
        if env_results_dir:
            return Path(env_results_dir)
        return base_dir / "results"

    def _resolve_auth_token(self, auth_token: str | None) -> str:
        """Resolve the auth token: explicit arg > generated token > AUTH_TOKEN env."""
        if auth_token:
            return auth_token
        token = get_token()
        if token:
            logger.info("Using generated auth token")
            return token
        token = environ.get("AUTH_TOKEN")
        if token:
            logger.info("Using auth token from environment variables")
            return token
        msg = "Failed to get authentication token"
        raise ValueError(msg)

    async def _process_samples(
        self,
        samples: list[tuple[str, str]],
        results_dir: Path,
        report_dir: Path,
        auth_token: str,
        existing_pdfs: set[str],
    ) -> tuple[int, int, int]:
        """Process all samples; return (success, skip, fail) counts."""
        success_count = 0
        skip_count = 0
        fail_count = 0

        async with aiohttp.ClientSession() as session:
            for lid, barcode in samples:
                try:
                    logger.info(f"Processing sample {lid} with barcode {barcode}")

                    if lid in existing_pdfs:
                        logger.info(f"PDF already exists for {lid}, skipping")
                        skip_count += 1
                        continue

                    result = await self.process_sample(
                        session,
                        lid,
                        barcode,
                        results_dir,
                        report_dir,
                        auth_token,
                        existing_pdfs,
                    )
                    if result:
                        logger.info(f"Successfully processed sample {lid}")
                        success_count += 1
                    else:
                        fail_count += 1
                    # Add a small delay between samples to avoid overwhelming the API
                    await asyncio.sleep(1)
                except aiohttp.ClientResponseError as e:
                    if e.status == _HTTP_UNAUTHORIZED:
                        logger.error("Authentication failed. Please check your auth token. Stopping processing.")
                        break
                    logger.error(f"Failed to process sample {lid}: {e}")
                    fail_count += 1
                except Exception as e:  # noqa: BLE001
                    logger.error(f"Failed to process sample {lid}: {e}")
                    fail_count += 1

        return success_count, skip_count, fail_count

    @staticmethod
    def _log_summary(total_samples: int, success: int, skip: int, fail: int, total_pdfs: int) -> None:
        """Log the PDF generation summary."""
        logger.info("PDF Generation Summary:")
        logger.info(f"- Total samples: {total_samples}")
        logger.info(f"- Generated PDFs in this run: {success}")
        logger.info(f"- Skipped samples (PDFs already exist): {skip}")
        logger.info(f"- Failed samples: {fail}")
        logger.info(f"- Total PDFs in output directory: {total_pdfs}")

    async def generate_reports(
        self,
        batch_id: str,
        data_dir: str | None,
        results_dir: str | None,
        base_dir: str | Path,
        report_dir: str | Path,
        auth_token: str | None = None,
    ) -> None:
        """Main function to extract sample info and generate PDF reports."""
        try:
            # Load environment variables from .env file
            dotenv.load_dotenv()

            # Convert string paths to Path objects
            base_dir = Path(base_dir)
            report_dir = Path(report_dir)

            # Setup metadata path — used to locate Excel/TSV sample manifest
            # Priority: explicit data_dir arg > DATA_DIR env var > base_dir/data
            data_dir_path = self._resolve_data_dir(data_dir, base_dir)
            metadata_dir = data_dir_path / "metadata"
            tsv_path = metadata_dir / f"{batch_id}.tsv"

            # Setup results path — used to locate sample JSON files for PDF generation
            # Priority: explicit results_dir arg > RESULTS_DIR env var > base_dir/results
            results_dir_path = self._resolve_results_dir(results_dir, base_dir)

            # Get authentication token
            auth_token = self._resolve_auth_token(auth_token)
            self.auth_token = auth_token

            # Find data file (Excel or TSV)
            data_file = self.find_data_file(batch_id, metadata_dir)
            if not data_file:
                msg = f"No Excel or TSV file found for batch {batch_id} in {metadata_dir}"
                raise FileNotFoundError(msg)  # noqa: TRY301

            logger.info(f"Found data file: {data_file}")

            # Create output directory
            report_dir.mkdir(parents=True, exist_ok=True)

            # Step 1: Extract sample info from data file
            logger.info(f"Extracting sample information from {data_file}")
            samples = self.extract_sample_info(str(data_file), str(tsv_path))

            # Get existing PDFs to avoid reprocessing
            existing_pdfs = self.get_existing_pdfs(report_dir)
            if existing_pdfs:
                logger.info(f"Found {len(existing_pdfs)} existing PDFs in {report_dir}")

            # Step 2: Process samples to generate PDFs
            success_count, skip_count, fail_count = await self._process_samples(
                samples,
                results_dir_path,
                report_dir,
                auth_token,
                existing_pdfs,
            )

            # Count total PDFs in output directory
            total_pdfs = len(list(report_dir.glob("*.pdf")))
            self._log_summary(len(samples), success_count, skip_count, fail_count, total_pdfs)

        except Exception as e:
            logger.error(f"Error generating reports: {e}")
            raise


def main() -> None:
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Generate PDF reports for mtDNA samples")
    parser.add_argument(
        "-b",
        "--batch_id",
        required=True,
        help="Batch ID (e.g., 20250210_mtDNA_17)",
    )
    parser.add_argument(
        "-d",
        "--data_dir",
        help="Path to data directory containing metadata/ (for sample manifest lookup). "
        "Falls back to DATA_DIR env var or base_dir/data if not provided.",
    )
    parser.add_argument(
        "-r",
        "--results_dir",
        required=True,
        help="Path to results directory containing sample JSON files for PDF generation "
        "(e.g., results/manual_pipeline/regenerate/{BATCH_ID}/)",
    )
    parser.add_argument(
        "--report_dir",
        required=True,
        help="Path to output directory for PDF reports",
    )
    parser.add_argument(
        "--auth_token",
        help="Authentication token for the API (if not provided, will try to get from .env or generate new one)",
    )

    args = parser.parse_args()

    # Use base directory to construct paths
    base_dir = Path(__file__).resolve().parent.parent

    # Create ReportGenerator instance and run
    generator = ReportGenerator()
    asyncio.run(
        generator.generate_reports(
            args.batch_id,
            args.data_dir,
            args.results_dir,
            base_dir,
            args.report_dir,
            args.auth_token,
        ),
    )


if __name__ == "__main__":
    main()
