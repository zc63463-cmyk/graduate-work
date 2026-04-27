#!/usr/bin/env python3
"""
提取两篇PDF论文的文本内容
"""
import sys
import os

# 添加PDF处理库
try:
    from pypdf import PdfReader
    print("Using pypdf library")
except ImportError:
    try:
        from PyPDF2 import PdfReader
        print("Using PyPDF2 library")
    except ImportError:
        print("Error: Please install pypdf or PyPDF2")
        print("Run: pip install pypdf")
        sys.exit(1)

def extract_pdf_text(pdf_path, output_path, max_pages=None):
    """提取PDF文本内容并保存到文件"""
    try:
        reader = PdfReader(pdf_path)
        num_pages = len(reader.pages)
        
        print(f"\nProcessing: {os.path.basename(pdf_path)}")
        print(f"Total pages: {num_pages}")
        
        if max_pages:
            num_pages = min(num_pages, max_pages)
            print(f"Extracting first {num_pages} pages...")
        
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(f"# PDF Extraction: {os.path.basename(pdf_path)}\n\n")
            f.write(f"Total pages: {len(reader.pages)}\n\n")
            f.write("---\n\n")
            
            for i in range(num_pages):
                page = reader.pages[i]
                text = page.extract_text()
                if text:
                    f.write(f"## Page {i+1}\n\n")
                    f.write(text)
                    f.write("\n\n")
        
        print(f"Text extracted to: {output_path}")
        return True
        
    except Exception as e:
        print(f"Error processing {pdf_path}: {e}")
        return False

def main():
    base_dir = r"C:\Users\20564\Desktop\Graduate\论文收集"
    
    papers = [
        {
            "pdf": os.path.join(base_dir, "源文件", "355853.355859.pdf"),
            "output": os.path.join(base_dir, "paper1_extracted.txt")
        },
        {
            "pdf": os.path.join(base_dir, "源文件", "355853.355865.pdf"),
            "output": os.path.join(base_dir, "paper2_extracted.txt")
        }
    ]
    
    for paper in papers:
        if os.path.exists(paper["pdf"]):
            extract_pdf_text(paper["pdf"], paper["output"])
        else:
            print(f"File not found: {paper['pdf']}")

if __name__ == "__main__":
    main()
