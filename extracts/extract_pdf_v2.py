from pypdf import PdfReader
import pytesseract
from pdf2image import convert_from_path
from PIL import Image
import os

pdf_path = r"C:\Users\20564\Desktop\Graduate\论文收集\直接文本补充\A MODIFIED GMRES METHOD FO...ONSYMMETRIC LINEAR SYSTEMS_牛强.pdf"
output_txt = r"C:\Users\20564\Desktop\Graduate\论文收集\直接文本补充\GMRES_paper_text.txt"

# 首先尝试用pypdf提取
print("尝试使用pypdf提取文本...")
try:
    reader = PdfReader(pdf_path)
    full_text = ""
    
    for i, page in enumerate(reader.pages):
        text = page.extract_text()
        full_text += f"\n\n=== 第 {i+1} 页 ===\n\n"
        full_text += text if text else f"[第{i+1}页无法提取文本，可能需要OCR]"
    
    # 检查提取的质量
    if "?" in full_text or len(full_text.strip()) < 100:
        print("pypdf提取效果不好，尝试OCR...")
        
        # 转换PDF为图像并使用OCR
        print("正在将PDF转换为图像（这可能需要一些时间）...")
        images = convert_from_path(pdf_path, dpi=300)
        
        full_text = ""
        for i, img in enumerate(images):
            print(f"正在OCR处理第 {i+1}/{len(images)} 页...")
            text = pytesseract.image_to_string(img, lang='eng+chi_sim')
            full_text += f"\n\n=== 第 {i+1} 页 ===\n\n"
            full_text += text
        
except Exception as e:
    print(f"错误: {e}")
    print("尝试直接使用OCR...")

# 保存提取的文本
with open(output_txt, 'w', encoding='utf-8') as f:
    f.write(full_text)

print(f"\n文本已保存到: {output_txt}")
print(f"\n前1000字符预览:")
print(full_text[:1000])
