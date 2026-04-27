import sys
import io
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8')

from pypdf import PdfReader

pdf_path = 'c:/Users/20564/Desktop/Graduate/论文收集/GMRES_translated.pdf'
reader = PdfReader(pdf_path)

full_text = ""
for i, page in enumerate(reader.pages):
    text = page.extract_text()
    if text:
        full_text += f"\n\n{'='*80}\n第 {i+1} 页\n{'='*80}\n\n"
        full_text += text

# 保存到文件
output_path = 'c:/Users/20564/Desktop/Graduate/论文收集/GMRES_full_text.txt'
with open(output_path, 'w', encoding='utf-8') as f:
    f.write(full_text)

print(f"已提取 {len(reader.pages)} 页内容")
print(f"保存到: {output_path}")
print(f"\n前3000字符预览:\n")
print(full_text[:3000])
