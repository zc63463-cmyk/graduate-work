import pdfplumber
import sys

pdf_path = r"C:\Users\20564\Desktop\Graduate\论文收集\直接文本补充\A MODIFIED GMRES METHOD FO...ONSYMMETRIC LINEAR SYSTEMS_牛强.pdf"

try:
    with pdfplumber.open(pdf_path) as pdf:
        print(f"总页数: {len(pdf.pages)}\n")
        print("=" * 80)
        
        full_text = ""
        for i, page in enumerate(pdf.pages):
            text = page.extract_text()
            if text:
                full_text += f"\n\n=== 第 {i+1} 页 ===\n\n"
                full_text += text
        
        # 保存到文件
        output_path = r"C:\Users\20564\Desktop\Graduate\论文收集\直接文本补充\GMRES_paper_extracted.txt"
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(full_text)
        
        print(f"PDF内容已提取到: {output_path}")
        print("\n前2000字符预览:")
        print(full_text[:2000])
        
except Exception as e:
    print(f"错误: {e}")
    sys.exit(1)
