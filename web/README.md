# Web UI cho mtDNA pipeline

Thay cho việc sửa tay `.env` và `scripts/batch_pipeline.sh` trước mỗi lần chạy.
Chọn quy trình trên trình duyệt, xem trước đúng dòng lệnh sẽ chạy, bấm Chạy, xem log trực tiếp.

## Chạy

```bash
cd /home/minhtq/mtdna-ui
bash web/start.sh              # cổng 8765, mở trong LAN
bash web/start.sh --port 9000
```

Không cần cài thêm gói nào — chỉ dùng thư viện chuẩn của Python 3.
Cố ý không dùng FastAPI: `.venv` của repo là venv đang chạy production, không nên cài thêm vào đó.

## Trang có gì

| Phần | Việc |
|---|---|
| Quy trình có sẵn | 5 nút điền sẵn theo SOP. Điền xong vẫn sửa lại từng ô được. |
| Các bước chạy | Tick tự do bước 1–6, không bị buộc theo quy trình. |
| Batch | Gõ mã batch, mỗi dòng một mã. `MS_110826_015` hoặc `20241223_mtDNA_12`. |
| Cấu hình | Engine, chế độ rerun, thư mục Sequencher / kết quả, chọn sheet, bật-tắt upload. |
| Xem trước | Thanh dưới cùng: thư mục thực tế, sheet, dòng lệnh chính xác, ZIP nào sẽ dùng, bao nhiêu file AB1, cảnh báo. |

Xem trước tự cập nhật mỗi khi đổi ô nào đó. Nếu nhập sai, nút Chạy bị khoá kèm lý do.

## Trong lúc chạy

Bấm Chạy thì form nhường chỗ cho **bảng tiến độ**: mỗi bước một dòng, kèm batch
đang xử lý.

```
✓ Chuẩn bị dữ liệu AB1        MS_270726_003
◐ Chạy engine phân tích       MS_270726_003
○ So sánh với Sequencher      MS_270726_003
○ Gộp kết quả so sánh
```

`✓` xong · `◐` đang chạy · `○` chưa tới · `—` bỏ qua · `✕` hỏng.

Log không đổ ra màn hình. Cảnh báo và lỗi được tách riêng thành một mục, vì đó
là phần thật sự cần đọc. Toàn bộ log vẫn còn, nằm trong mục "Xem log chi tiết"
thu gọn, và ghi ra `web/runs/<id>.log`.

Bảng tiến độ dựng ở server từ chính dòng log pipeline vẫn in — không phải thêm
thiết bị đo nào vào `pipeline.sh`. Nếu sau này log đổi chữ, sửa các mẫu nhận
dạng trong lớp `Progress` ở `web/server.py`.

## Google Sheet

Mỗi nhiệm vụ một sheet riêng, khai trong `.env`:

| Biến | Dùng cho |
|---|---|
| `GOOGLE_SHEET_ID` | batch thường quy |
| `GOOGLE_SHEET_ID_TRAINING` | batch training |
| `GOOGLE_SHEET_ID_RERUN` | batch rerun |

Trên UI chọn theo tên. **ID không hiện ở bất cứ đâu trên trang** — không trong ô
chọn, không trong dòng lệnh xem trước, không trong log. Nó nằm trong `.env` và đi
thẳng từ server tới tiến trình con; dòng lệnh hiển thị là `--sheet-id ***`.

Cả file `web/runs/*.log` cũng đã che, vì file log sống lâu hơn phiên làm việc.

Phần xem trước vẫn ghi rõ **tên** sheet sẽ ghi vào trước khi chạy.

Muốn dùng một sheet khác thì thêm biến vào `.env` rồi khai trong `SHEET_KEYS`
ở `web/server.py`; dòng lệnh vẫn nhận `--sheet-id` như cũ.

**Bật upload thì bắt buộc phải chọn sheet** — không có chuyện bỏ trống rồi rơi về
sheet chính. Đây đúng là kiểu tai nạn đã từng xảy ra: một lượt rerun ghi đè
worksheet của batch thường quy.

Sheet nào chưa khai trong `.env` thì không hiện trên UI.

### Worksheet nào được tạo

| Chế độ | Worksheet |
|---|---|
| Thường quy | `Batch_<mã batch>` cho từng batch, cộng thêm `All_Batches` gộp lại |
| Rerun | `All_Batches`, **cộng thêm** một tab mang tên thư mục chứa ZIP, ví dụ `Rerun_260801_RT_013`. Không tạo tab cho từng batch. |

Rerun gom vài mẫu từ rất nhiều batch, nên mỗi batch một tab sẽ chôn cả sheet
dưới đống tab gần như rỗng. Tab mang tên thư mục nguồn giữ lại kết quả của từng
đợt rerun — `All_Batches` thì luôn bị đợt sau ghi đè.

### Định dạng worksheet

Uploader trước đây chỉ ghi giá trị. `worksheet.clear()` xoá dữ liệu chứ **không**
xoá định dạng, nên worksheet có sẵn giữ nguyên format cũ còn worksheet mới tạo
thì trắng trơn theo mặc định Google — cùng một pipeline mà hai sheet trông khác
hẳn nhau.

Giờ mỗi lần upload đều áp một bộ định dạng chung:

- đóng băng hàng tiêu đề
- tiêu đề in đậm, nền `D9EAF7`, canh giữa — cùng màu với cách phần còn lại của
  dự án định dạng workbook (`autosize_and_style()` trong
  `src/tools/mutation_surveyor/merge_final_profiles_vs_truth.py`)
- bật bộ lọc
- **mọi cột một bề rộng: 90px**, ô bật xuống dòng

Bề rộng đồng đều là cố ý. Co từng cột vừa nội dung thì bảng lởm chởm — một cột
rộng gấp sáu lần cột bên cạnh, mắt lạc dòng. Cố định bề rộng rồi cho chữ chảy
xuống thì bảng thành lưới đều, và 13 cột gói trong 1170px, vừa một màn hình.

Đổi bằng `--column-width PX` hoặc sửa `UNIFORM_COLUMN_WIDTH_PX` trong
`src/modules/sheets/uploader.py`. Đo trên batch 90 mẫu thật:

| Bề rộng | Ký tự/dòng | Tổng ngang | Dòng cao nhất | Dòng trung bình |
|---|---|---|---|---|
| **90px** | 11 | **1170px** | 594px | 189px |
| 120px | 16 | 1560px | 486px | 116px |
| 160px | 21 | 2080px | 306px | 84px |
| 200px | 27 | 2600px | 288px | 66px |

Hẹp thì gọn ngang nhưng cao dọc — không tránh được, chỉ chọn điểm cân bằng.

Chiều cao dòng trả về cho Google tự co (`autoResizeDimensions`) — dòng giữ
nguyên chiều cao đã đặt trước đó, không tự co lại, nên phải nói rõ.

Tất cả gộp trong **một** lần gọi API.

Định dạng hỏng thì chỉ ghi cảnh báo, không làm hỏng lần upload — dữ liệu quan
trọng hơn.

Muốn giữ nguyên format cũ của một worksheet thì thêm `--no-format` khi gọi
`uploader.py`. Lưu ý: chạy có định dạng sẽ **ghi đè** các chỉnh tay trước đó
(bề rộng cột, màu) trên worksheet đó.

Mặc định của từng quy trình:

| Quy trình | Sheet | Upload |
|---|---|---|
| Automate + compare, Compare lại | chính | bật |
| Training | training | bật |
| Rerun | rerun | bật |
| Xuất FASTA | không dùng | tắt |

Bật hay tắt vẫn đổi được từng lần chạy trên UI, hoặc bằng `--upload` /
`--no-upload` ở dòng lệnh.

## Chạy nhiều batch cùng lúc

Một đợt rerun thường gom vài mẫu từ rất nhiều batch, chạy lần lượt thì lâu.

**Mặc định là tự đếm**: nhập bao nhiêu batch thì chạy bấy nhiêu cùng lúc. Không
có lý do giữ batch lại chờ, và trần CPU đã chặn tổng khối lượng rồi. Muốn ép
xuống thì điền số vào ô **Số batch chạy cùng lúc**, hoặc `-j N` ở dòng lệnh:

```bash
bash scripts/batch_pipeline.sh --rerun -s 3,4 -p blastn -j 4 -b BATCH1,BATCH2,...
```

Đo trên 4 batch training, bước 3+4, engine tracy:

| | Thời gian |
|---|---|
| lần lượt | 189s |
| 4 batch cùng lúc | **55s** |

Kết quả giống hệt nhau — 8/8 file so sánh trùng từng byte.

An toàn vì mọi đường dẫn của một batch đều mang mã batch của nó
(`data/raw/<batch>`, `results/tools/*/<batch>`, `modules/comparison/<batch>`).
Phần dùng chung là bước gộp, và nó chạy sau khi mọi batch xong.

### Phân bổ CPU

Để mặc thì mỗi công cụ lấy `os.cpu_count()` — cả 256 nhân — nên một lượt chạy
chiếm trọn máy và N batch song song chỉ giành nhau. Mỗi lượt giờ được cấp
**128 nhân**, chia đều cho các batch chạy cùng lúc trong lượt đó:

| Batch cùng lúc | Worker mỗi batch |
|---|---|
| 1 | 128 |
| 4 | 32 |
| 8 | 16 |
| 23 | 5 |

Xem thêm [Chia CPU giữa các lượt](#chia-cpu-giữa-các-lượt) về phần nhiều người
chạy cùng lúc.

Log của mỗi batch được gắn tiền tố mã batch, nếu không thì các luồng trộn vào
nhau và không biết dòng `[Step 3/6]` là của batch nào.

## Sửa danh sách quy trình

`web/presets.json`. Sửa xong tải lại trang là có hiệu lực, không cần khởi động lại server.

## Vài điểm về cách hoạt động

- Lệnh chạy luôn là mảng tham số, **không bao giờ qua shell** — không thể chèn lệnh qua ô nhập.
- Mã batch, tên bước, engine, sheet ID đều bị kiểm định dạng trước khi chạy.
- Nhiều người chạy cùng lúc được — xem mục dưới.
- Nút Dừng gửi SIGTERM cho cả cây tiến trình con.
- Log của từng lần chạy được ghi lại ở `web/runs/<id>.log`.
- Tải lại trang giữa chừng vẫn bám lại được lần chạy đang chạy, bảng tiến độ vẽ lại đúng chỗ đang đứng.
- Không có đăng nhập. Ai vào được LAN là chạy được pipeline — đừng mở ra Internet.

## Nhiều người dùng cùng lúc

Không giới hạn số lượt chạy song song. Chặn duy nhất là **cùng một batch**: hai
lượt chạy chung một batch sẽ xoá output của nhau giữa chừng, nên lượt sau bị từ
chối kèm tên người đang chạy nó.

```
Batch MS_270726_003 đang được chạy ở một lượt khác (Minh).
Đợi lượt đó xong, hoặc bỏ batch đó ra.
```

Ô **Người chạy** trên form chỉ để biết lượt nào của ai — không có đăng nhập,
không phân quyền. Tên được nhớ trong trình duyệt, khỏi gõ lại. Nút **Các lượt
đang chạy** ở đầu trang liệt kê tất cả, kèm nút Xem và Dừng cho từng lượt. Dừng
lượt của người khác thì có hỏi lại, nêu rõ tên và batch.

Tải lại trang giữa chừng: một lượt đang chạy thì bám thẳng vào nó, nhiều lượt
thì hiện danh sách để chọn.

Đo thật, hai người hai lượt cùng lúc, mỗi lượt 2 batch: 51s và 53s, cả hai đều
`rc=0`. Chạy chung một thư mục kết quả cũng được — bốn thư mục comparison còn
nguyên vẹn, riêng `all_batches_comparison` là của người xong sau.

### Chia CPU giữa các lượt

**Mỗi lượt 128 nhân, ai cũng như nhau.** Trên máy 256 nhân, ba người chạy cùng
lúc sẽ xin tổng 384 — vượt, và điều đó là cố ý.

Từng cân nhắc chia khẩu phần theo ngân sách chung, nhưng số worker của một lượt
bị ấn định lúc khởi động và không co lại được cho người tới sau. Chia khẩu phần
tức là ai vào sau thì chạy chậm suốt thời gian lượt trước còn chạy — cái giá đó
đắt hơn việc để hệ điều hành tự xoay khi bị đăng ký vượt.

`/api/status` vẫn báo `used`, `total` và cờ `overcommitted` để biết máy đang
căng tới đâu.

Đổi mức bằng `MTDNA_MAX_CORES` (mặc định 128). Chạy từ dòng lệnh không qua UI
thì dùng `--max-cores N`.

## Tự kiểm

```bash
make ui-test              # ~10 giây, chỉ kiểm API và giao diện
make ui-test FULL=1       # thêm 2 lượt pipeline thật, ghi vào thư mục tạm
```

Bản `FULL=1` chạy pipeline thật nhưng **không gửi gì lên Google**. Hai lớp chặn:

- Các lượt chạy pipeline đều bị ép `--no-upload`, và test tự dừng nếu dựng ra
  một dòng lệnh có khả năng chạm tới sheet.
- Riêng phần kiểm tra "upload có gọi đúng sheet không" thì chặn ngay ở
  `uploader.py` bằng một `python` giả đặt trước trong PATH: nó ghi lại lời gọi
  rồi trả về 0, mọi lệnh python khác chuyển tiếp cho python thật. Nhờ vậy vẫn
  khẳng định được nó *định* ghi vào sheet nào mà không thật sự ghi.

Kết quả ghi vào thư mục tạm, không đụng `results` trên NAS.

## Quan hệ với các file khác

UI chỉ gọi `scripts/batch_pipeline.sh` với các cờ dòng lệnh, không ghi đè `.env`.
Giá trị nào không điền trên form thì `batch_pipeline.sh` lấy từ `.env` như cũ.
Các script `run*.sh` ở gốc repo vẫn dùng được song song, làm cùng một việc từ terminal.
