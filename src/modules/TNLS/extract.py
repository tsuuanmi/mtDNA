#!/usr/bin/env python
import re
from os import environ
from pathlib import Path

import dotenv
import pandas as pd
from loguru import logger

# Load environment variables from .env file
dotenv.load_dotenv()

# Use DATA_DIR and RESULTS_DIR from .env if they exist, otherwise use default paths
DATA_DIR = Path(environ.get("DATA_DIR"))  # type: ignore[reportArgumentType]
RESULTS_DIR = Path(environ.get("RESULTS_DIR"))  # type: ignore[reportArgumentType]



class SoldierNameExtractor:
    """Class for extracting soldier names from notes in a dataframe."""

    def __init__(self) -> None:
        # Different patterns to match soldier names
        self.patterns = [
            # Standard patterns for direct mentions
            r"(?:liệt sĩ|Liệt sĩ|LS|liệt sỹ|Liệt sỹ)\s+([\w\s,]+?)(?:$|\.|\s+và|\s+con|/)",
            # Pattern for 'của liệt sĩ' format - expanded with "sỹ" variant
            r"của\s+(?:liệt sĩ|Liệt sĩ|LS|liệt sỹ|Liệt sỹ|chị\s+LS)\s+([\w\s,]+?)(?:$|\.|\s+và|\s+\()",
            # Pattern for 'bên họ' format
            r"bên họ\s+(?:liệt sĩ|Liệt sĩ|LS|liệt sỹ|Liệt sỹ)\s+([\w\s,]+?)(?:$|\.|\s+và)",
            # Simplified pattern for just the identifier, with commas
            r"(?:liệt sĩ|Liệt sĩ|LS|liệt sỹ|Liệt sỹ)\s+([\w\s,]+)",
            # Pattern for "Con của chị LS" and similar formats
            r"(?:Con|Em|Cháu)\s+(?:của|gái|trai|ruột)?\s+(?:chị|anh|em)?\s+(?:LS|liệt sĩ|liệt sỹ)\s+([\w\s,]+?)(?:$|\.|\s+và|\s+\()",  # noqa: E501
            # Handle colon (:) after relationship and LS - expanded with more relationship terms
            r"(?:Con|Em|Cháu|Mẹ đẻ|Anh|Chị|Mẹ|Ba|Bố|Cha|Cháu nội|Cháu ngoại|Anh trai|Chị gái|Em gái|Em trai|Chú|Bác|Thím|Dì)\s+(?:gái|trai|ruột|đẻ|họ|nội|ngoại)?\s+(?:LS|liệt sĩ|liệt sỹ):\s+([\w\s,()]+?)(?:$|\.|\s+và)",  # noqa: E501
            # Format with space before colon - "Anh trai LS : Nguyễn Văn Nhỏ"
            r"(?:Con|Em|Cháu|Mẹ đẻ|Anh|Chị|Mẹ|Ba|Bố|Cha|Cháu nội|Cháu ngoại|Anh trai|Chị gái|Em gái|Em trai|Chú|Bác|Thím|Dì)\s+(?:gái|trai|ruột|đẻ|họ|nội|ngoại)?\s+(?:LS|liệt sĩ|liệt sỹ)\s+:\s+([\w\s,()]+?)(?:$|\.|\s+và)",  # noqa: E501
            # Format like "Em LS: Bùi Văn Hưng" - expanded with more relationship terms
            r"(?:Con|Em|Cháu|Mẹ đẻ|Anh|Chị|Mẹ|Ba|Bố|Cha|Cháu nội|Cháu ngoại|Anh trai|Chị gái|Em gái|Em trai|Chú|Bác|Thím|Dì)\s+(?:LS|liệt sĩ|liệt sỹ):\s+([\w\s,()]+?)(?:$|\.|\s+và)",  # noqa: E501
            # Format like "Em Ls Phan Huy Sơn" (lowercase 'Ls') - FIXED: Ensure no space requirement
            r"(?:Con|Em|Cháu|Mẹ đẻ|Anh|Chị|Mẹ|Ba|Bố|Cha|Cháu nội|Cháu ngoại|Anh trai|Chị gái|Em gái|Em trai|Chú|Bác|Thím|Dì)\s+(?:Ls|ls|SL)\s*([\w\s,()]+?)(?:$|\.|\s+và)",  # noqa: E501
            # Format like "Em ruột Ls: Hà Bình Minh"
            r"(?:Con|Em|Cháu|Mẹ đẻ|Anh|Chị|Mẹ|Ba|Bố|Cha|Cháu nội|Cháu ngoại|Anh trai|Chị gái|Em gái|Em trai|Chú|Bác|Thím|Dì)\s+(?:gái|trai|ruột|đẻ|họ|nội|ngoại)?\s+(?:Ls|ls):\s+([\w\s,()]+?)(?:$|\.|\s+và)",  # noqa: E501
            # Format like "em gái ruột LS: Nguyễn Trọng Mậu" (lowercase first word)
            r"(?:con|em|cháu|mẹ đẻ|anh|chị|mẹ|ba|bố|cha|cháu nội|cháu ngoại|anh trai|chị gái|em gái|em trai|chú|bác|thím|dì)\s+(?:gái|trai|ruột|đẻ|họ|nội|ngoại)?\s+(?:LS|liệt sĩ|liệt sỹ):\s+([\w\s,()]+?)(?:$|\.|\s+và)",  # noqa: E501
            # Format for "Em SL Trần Xuân Trình" (SL instead of LS) - FIXED: Made asterisk after space
            r"(?:Con|Em|Cháu|Mẹ đẻ|Anh|Chị|Mẹ|Ba|Bố|Cha|Cháu nội|Cháu ngoại|Anh trai|Chị gái|Em gái|Em trai|Chú|Bác|Thím|Dì)\s+(?:gái|trai|ruột|đẻ|họ|nội|ngoại)?\s*(?:SL)\s*([\w\s,()]+?)(?:$|\.|\s+và)",  # noqa: E501
            # "Thân nhân liệt sữ" variant (typo)
            r"(?:Thân nhân|Người thân)\s+(?:liệt sữ|liệt sỹ|Liệt sữ|Liệt sỹ)\s+([\w\s,]+?)(?:$|\.|\s+và)",
            # Pattern for "Cháu ruột Nguyễn Hữu Tấc, con em trai liệt sĩ"
            r"(?:Con|Em|Cháu|Mẹ đẻ|Anh|Chị|Mẹ|Ba|Bố|Cha|Cháu nội|Cháu ngoại|Anh trai|Chị gái|Em gái|Em trai|Chú|Bác|Thím|Dì)\s+(?:gái|trai|ruột|đẻ|họ|nội|ngoại)?\s+([\w\s,()]+),?\s+(?:con|em).*?(?:liệt sĩ|LS|liệt sỹ)",  # noqa: E501
            # Format like "Con của chị gái Nguyễn Dũng Đàng" - relationship without LS
            r"(?:Con|Cháu)\s+của\s+(?:chị|anh|em)\s+(?:gái|trai|ruột)?\s+([\w\s,]+)(?:$|\.)",
            # Format like "Cháu (con của chị gái) Nguyễn Khắc Mễ"
            r"(?:Con|Em|Cháu|Mẹ đẻ|Anh|Chị|Mẹ|Ba|Bố|Cha)\s+\((?:con|cháu|em|anh|chị)\s+của\s+(?:chị|anh|em|bố|mẹ)\s+(?:gái|trai|ruột)?\)\s+([\w\s,]+)(?:$|\.)",  # noqa: E501
            # Handle cases without LS indicator (simplified)
            r"^(?:Con|Em|Cháu|Mẹ đẻ|Anh|Chị|Mẹ|Ba|Bố|Cha|Cháu nội|Cháu ngoại|Anh trai|Chị gái|Em gái|Em trai|Chú|Bác|Thím|Dì)\s+(?:gái|trai|ruột|đẻ|họ|nội|ngoại)\s+([\w\s,()]+)$",  # noqa: E501
            # "Em và chị LS: Nguyễn Văn Mão, Nguyễn Văn Minh"  # noqa: ERA001
            r"(?:Con|Em|Cháu|Anh|Chị)\s+và\s+(?:con|em|cháu|anh|chị)\s+(?:LS|liệt sĩ|liệt sỹ|Ls|ls):\s+([\w\s,()]+?)(?:$|\.)",  # noqa: E501
            # "Chị/Em LS (cùng mẹ) Phạm Văn Hiến" - FIXED: Made space optional
            r"(?:Con|Em|Cháu|Anh|Chị)/(?:con|em|cháu|anh|chị)\s+(?:LS|liệt sĩ|liệt sỹ|Ls|ls)(?:\s+\([^)]*\))?\s*([\w\s,()]+?)(?:$|\.)",  # noqa: E501
            # "Em gái của LS: Nguyễn Minh Trí"  # noqa: ERA001
            r"(?:Con|Em|Cháu|Anh|Chị)\s+(?:gái|trai|ruột)?\s+của\s+(?:LS|liệt sĩ|liệt sỹ|Ls|ls):\s+([\w\s,()]+?)(?:$|\.)",  # noqa: E501
            # "Con chị gái Vũ Ngọc Đáng" (without "của")
            r"(?:Con|Cháu)\s+(?:chị|anh|em)\s+(?:gái|trai|ruột)?\s+([\w\s,]+)(?:$|\.)",
            # "Mẹ Ls Hồ Hữu Lạc" (without a space, directly attached) - FIXED: Changed to fix missing spaces
            r"(?:Mẹ|Ba|Bố|Cha|Anh|Chị|Em)(?:Ls|LS|ls|SL)\s*([\w\s,()]+?)(?:$|\.)",
            # "Em của LS: Đặng Đình Huyền"  # noqa: ERA001
            r"(?:Con|Em|Cháu|Anh|Chị|Mẹ|Ba|Bố|Cha)\s+của\s+(?:LS|liệt sĩ|liệt sỹ|Ls|ls):\s+([\w\s,()]+?)(?:$|\.)",
            # SPECIAL CASE: "Con gái ruột LS: Trần Cao Thanh" - FIXED: Added specific pattern
            r"Con\s+gái\s+ruột\s+LS:\s+([\w\s,()]+?)(?:$|\.)",
            # Special case for "Cha đẻ của liệt sĩ." - explicit extraction with empty name
            r"Cha\s+đẻ\s+của\s+liệt sĩ\.",
            # Pattern to handle multiple LS mentions in sequence like "LS Name1, LS Name2"
            r"(?:liệt sĩ|Liệt sĩ|LS|liệt sỹ|Liệt sỹ)\s+([\w\s,]+?)(?:,\s*(?:liệt sĩ|Liệt sĩ|LS|liệt sỹ|Liệt sỹ))",
            # Pattern to handle "Relationship LS Name1, LS Name2" format
            r"(?:Con|Em|Cháu|Mẹ đẻ|Anh|Chị|Mẹ|Ba|Bố|Cha|Cháu nội|Cháu ngoại|Anh trai|Chị gái|Em gái|Em trai|Chú|Bác|Thím|Dì)\s+(?:gái|trai|ruột|đẻ|họ|nội|ngoại)?\s+(?:LS|liệt sĩ|liệt sỹ|Ls|ls)\s+([\w\s\d]+?)(?:,\s*(?:LS|liệt sĩ|liệt sỹ|Ls|ls))",  # noqa: E501
            # Pattern to handle "Mẹ đẻ LS Name YYYY" format specifically
            r"(?:Mẹ đẻ|Mẹ|Ba|Bố|Cha)\s+(?:LS|liệt sĩ|liệt sỹ|Ls|ls)\s+([\w\s\d]+?)(?:$|\.|\s+và)",
            # Pattern to capture names with years (e.g., "Nguyễn Văn Chung 1958")
            r"(?:LS|liệt sĩ|liệt sỹ|Ls|ls)\s+([\w\s]+\s+\d{4})(?:$|\.|\s+và)",
        ]

        # Define patterns that indicate no useful information for extraction
        self.no_info_patterns = [
            "Không có thông tin",
            "Thiếu TTLS",
            "Thiếu TT LS",
            "TT về MQH chưa rõ",
            "CT 1000 mẫu miễn phí",
            "Mẫu tặng",
            "email sai, SDT sai",
            "Không có tt LS",
            "Không có tt",
            "Không có tên LS",
            "thông báo",
            "ngõ",
            "ngách",
            "Phường",
            "Quận",
            "Anh Thành xin duyệt chạy mẫu - cán bộ BCA",
            "không có thông tin LS",
            "KQ gửi về",
            "Mẫu gửi đi",
            "-",
            "Mẫu gửi đi",
            "Không rõ thông tin LS",
            "CC email sai",
            "Lỗi mail",
            "Không ghi",
            "Hội Thảo",
            "Lấy mẫu ở",
            "Mẫu KH liên hệ",
            "Làm thêm Diamond",
            "Làm thêm genemap",
            "Lấy mẫu tại",
            "KH nhận kết quả",
            "không xác thực được",
            "gia đình khiếm thính",
        ]

    def _clean_name(self, name: str) -> list[str]:
        """Clean a single soldier name by removing unwanted characters."""
        # Extract content within parentheses if it's potentially a name
        parentheses_content = re.findall(r"\(([\w\s]+)\)", name)
        additional_names = []

        # Check if content in parentheses is likely a name (not an ID)
        additional_names.extend(
            content for content in parentheses_content
            if not re.search(r"[A-Z0-9]{5,}", content) and not re.search(r"F\d+[A-Z]+", content)
        )

        # Remove parentheses content after extracting potential names
        name = re.sub(r"\s*\([^)]*\)", "", name)

        # Remove trailing slash
        name = re.sub(r"\s*/$", "", name)

        # Preserve years in names (like "Name 1958")
        # But remove other non-name info after the year
        name = re.sub(r"(\s+\d{4}).*$", r"\1", name)

        # Trim whitespace
        name = name.strip()

        # Return the main name and any additional names found in parentheses
        if additional_names:
            return [name, *additional_names]
        return [name]

    def _extract_names_from_note(self, note: str | None) -> list[str]:  # noqa: C901,PLR0912
        """Extract all soldier names from a single note."""
        # Handle empty or no-name cases
        if pd.isna(note):  # type: ignore[reportGeneralTypeIssues]
            return []

        # After pd.isna check, note is a string
        note_text = str(note)

        # Skip notes that explicitly state there's no information
        if any(pattern in note_text for pattern in self.no_info_patterns):
            return []

        if "Không ghi tên liệt sĩ" in note_text:
            return ["Không có tên"]

        all_names = []

        # Special handling for multiple LS mentions in the same note
        # First, try to find all occurrences of "LS Name" patterns
        ls_pattern = r"(?:liệt sĩ|Liệt sĩ|LS|liệt sỹ|Liệt sỹ)\s+([\w\s\d]+?)(?:$|\.|\s+và|,|\s+con)"
        ls_matches = re.finditer(ls_pattern, note_text)
        for match in ls_matches:
            name = match.group(1).strip()
            if name and len(name) > 1:
                all_names.append(name)

        # If we found multiple LS names, return them
        if len(all_names) > 1:
            return all_names

        # Otherwise, proceed with the regular pattern matching
        all_names = []
        for pattern in self.patterns:
            matches = re.findall(pattern, note_text)
            if matches:
                for match in matches:
                    # Skip empty matches or matches without a capturing group
                    if not match:
                        continue

                    # First split by "và" if present
                    if " và " in match:
                        parts = [p.strip() for p in match.split(" và ")]
                        for part in parts:
                            # Then check for commas in each part
                            if "," in part:
                                all_names.extend([p.strip() for p in part.split(",")])
                            else:
                                all_names.append(part)
                    # If no "và", check for commas
                    elif "," in match:
                        all_names.extend([p.strip() for p in match.split(",")])
                    else:
                        all_names.append(match.strip())

                # If we found names with this pattern, no need to try others
                if all_names:
                    break

        # Special case handling for complex cases
        if not all_names and not any(pattern in note_text.lower() for pattern in self.no_info_patterns):
            # Try to extract the name that appears before "con" or "em" and after a relationship
            complex_pattern = r"(?:Cháu|Con|Em)\s+(?:ruột|gái|trai|nội|ngoại)?\s+([\w\s]+?)(?:,|\s+(?:con|em|là))"
            matches = re.findall(complex_pattern, note_text)
            if matches:
                all_names.extend([m.strip() for m in matches])

        # Clean up names and filter out very short ones
        clean_names = []
        for name in all_names:
            cleaned_results = self._clean_name(name)
            for cleaned in cleaned_results:
                if cleaned and len(cleaned) > 1 and cleaned not in clean_names:
                    clean_names.append(cleaned)

        return clean_names

    def extract_names(self, df: pd.DataFrame) -> pd.DataFrame:
        """Extract soldier names from the entire dataframe."""
        result_df = df.copy()

        # Initialize the Soldier_Name column
        result_df["Soldier_Name"] = None

        # Priority order: Liet si > Family > Note extraction

        # 1. First check if "Liet si" column exists (highest priority)
        if "Liet si" in result_df.columns:
            logger.info("Found 'Liet si' column, using as primary source")
            # Use Liet si column data where it exists and is not null/empty
            liet_si_mask = result_df["Liet si"].notna() & (result_df["Liet si"].str.strip() != "")
            result_df.loc[liet_si_mask, "Soldier_Name"] = result_df.loc[liet_si_mask, "Liet si"]
            logger.info(f"Used {liet_si_mask.sum()} soldier names from 'Liet si' column")

        # 2. Check if Family column exists and use where Liet si is empty
        if "Family" in result_df.columns:
            logger.info("Found 'Family' column, using where 'Liet si' is empty")
            # Use Family column data where Soldier_Name is still empty
            family_mask = (
                result_df["Family"].notna() & (result_df["Family"].str.strip() != "") & result_df["Soldier_Name"].isna()
            )
            result_df.loc[family_mask, "Soldier_Name"] = result_df.loc[family_mask, "Family"]
            logger.info(f"Used {family_mask.sum()} soldier names from Family column")

        # 3. Extract from Note column for remaining rows
        remaining_rows = result_df["Soldier_Name"].isna()
        if remaining_rows.any():  # type: ignore[reportGeneralTypeIssues]
            logger.info(f"Extracting soldier names from Note column for {remaining_rows.sum()} remaining rows")

            # Process remaining rows
            for i, row in result_df[remaining_rows].iterrows():
                names = self._extract_names_from_note(row["Note"])  # type: ignore[reportArgumentType]
                if names:
                    result_df.loc[i, "Soldier_Name"] = ", ".join(names)

        return result_df


def load_data(file_path: str | Path) -> pd.DataFrame:
    """Load data from a TSV file."""
    try:
        return pd.read_csv(file_path, sep="\t", low_memory=False)
    except (FileNotFoundError, pd.errors.ParserError, OSError) as e:
        logger.error(f"Error loading data from {file_path}: {e}")
        raise


def strip_decimal_zeros(df: pd.DataFrame) -> pd.DataFrame:
    """Strip '.0' from all columns in the dataframe."""
    df_cleaned = df.copy()

    for col in df_cleaned.columns:
        # Convert to string and strip .0 from numeric values
        df_cleaned[col] = df_cleaned[col].astype(str).str.replace(r"\.0$", "", regex=True)
        # Convert 'nan' back to empty string for better readability
        df_cleaned[col] = df_cleaned[col].replace("nan", "")

    return df_cleaned


def save_data(df: pd.DataFrame, output_path: str | Path) -> None:
    """Save dataframe to a TSV file."""
    try:
        # Strip .0 from all columns before saving
        df_cleaned = strip_decimal_zeros(df)

        # Exclude Family column from final output if it exists
        if "Family" in df_cleaned.columns:
            df_cleaned = df_cleaned.drop("Family", axis=1)
            logger.info("Excluded Family column from final output")

        df_cleaned.to_csv(output_path, sep="\t", index=False)
        logger.info(f"Data saved successfully to {output_path} (with .0 stripped from all columns)")
    except OSError as e:
        logger.error(f"Error saving data to {output_path}: {e}")
        raise


def evaluate_extraction(df: pd.DataFrame) -> tuple[int, int, float]:
    """Evaluate the success of name extraction."""
    success_count = df["Soldier_Name"].notna().sum()
    total_count = df.shape[0]
    success_rate = success_count / total_count if total_count > 0 else 0

    logger.info(f"Successfully extracted names: {success_count} out of {total_count} ({success_rate:.2%})")

    # Show breakdown of sources
    source_counts = {}

    # Count from Liet si column
    if "Liet si" in df.columns:
        liet_si_count = (df["Liet si"].notna() & (df["Liet si"].str.strip() != "")).sum()
        source_counts["Liet si"] = liet_si_count

    # Count from Family column
    if "Family" in df.columns:
        family_count = (df["Family"].notna() & (df["Family"].str.strip() != "")).sum()
        source_counts["Family"] = family_count

    # Count from Note extraction (remaining)
    note_extraction_count = success_count - sum(source_counts.values())
    source_counts["Note extraction"] = note_extraction_count

    logger.info(f"Source breakdown: {', '.join(f'{k}: {v}' for k, v in source_counts.items())}")

    # Check rows where extraction failed but should have succeeded
    # Get the no_info_patterns from an instance of SoldierNameExtractor
    extractor = SoldierNameExtractor()
    no_info_patterns = extractor.no_info_patterns

    # Create a filter to exclude notes that shouldn't be considered as failures
    # Use a different approach to avoid the regex warning
    no_info_mask = pd.Series([False] * len(df), index=df.index)
    for pattern in no_info_patterns:
        no_info_mask = no_info_mask | df["Note"].str.contains(pattern, case=False, na=False)

    no_info_filter = ~no_info_mask

    failed = df[
        (df["Soldier_Name"].isna())
        & (~df["Note"].str.contains("Không ghi tên liệt sĩ", na=True))
        & (df["Note"].notna())
        & no_info_filter
    ]

    logger.info(f"Failed extractions (excluding rows without soldier name and no-info notes): {failed.shape[0]}")

    # Display more examples of failures for debugging
    if not failed.empty:
        sample_failures = failed.sample(min(25, failed.shape[0]))
        logger.info("Sample failures:")
        for _, row in sample_failures.iterrows():
            logger.info(f"Note: '{row['Note']}'")

    return int(success_count), int(total_count), float(success_rate)  # type: ignore[reportCallIssue]


def main() -> None:

    input_file = DATA_DIR / "TNLS/STR_mtDNA_with_ids.tsv"
    output_file = DATA_DIR / "TNLS/TNLS_with_soldier_names.tsv"

    # Load data
    logger.info(f"Loading data from {input_file}")
    df = load_data(input_file)

    # Extract soldier names
    logger.info("Extracting soldier names")
    extractor = SoldierNameExtractor()
    df_with_names = extractor.extract_names(df)

    # Evaluate the extraction
    evaluate_extraction(df_with_names)

    # Save the results
    save_data(df_with_names, output_file)


if __name__ == "__main__":
    main()
