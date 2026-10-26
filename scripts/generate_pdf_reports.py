#!/usr/bin/env python3
"""
Generate PDF reports from markdown files using fpdf2.
Creates both Turkish and English versions.
"""

import os
import sys
from pathlib import Path
from fpdf import FPDF

class MarkdownToPDF(FPDF):
    def __init__(self):
        super().__init__()
        # Try to add DejaVu fonts
        font_paths = [
            '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf',
            '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf',
            'C:/Windows/Fonts/DejaVuSans.ttf',
            'C:/Windows/Fonts/DejaVuSans-Bold.ttf',
        ]
        for fp in font_paths:
            if os.path.exists(fp):
                if 'Bold' in fp:
                    self.add_font('DejaVu', 'B', fp, uni=True)
                else:
                    self.add_font('DejaVu', '', fp, uni=True)
        self.set_auto_page_break(auto=True, margin=20)
        
    def header(self):
        self.set_font('DejaVu', 'B', 10)
        self.cell(0, 5, 'Zamanda Yolculuk - AMR Brill-Lindquist Stability Report', 0, 1, 'C')
        self.line(10, self.get_y(), 200, self.get_y())
        self.ln(5)
        
    def footer(self):
        self.set_y(-15)
        self.set_font('DejaVu', '', 8)
        self.cell(0, 10, f'Sayfa {self.page_no()}/{{nb}}', 0, 0, 'C')

def read_markdown(filepath):
    """Read markdown file content."""
    with open(filepath, 'r', encoding='utf-8') as f:
        return f.read()

def markdown_to_pdf(md_content, output_path, title="Report"):
    """Convert markdown content to PDF."""
    pdf = MarkdownToPDF()
    pdf.alias_nb_pages()
    pdf.add_page()
    pdf.set_font('DejaVu', '', 10)
    
    lines = md_content.split('\n')
    in_code_block = False
    code_lines = []
    
    for line in lines:
        # Handle code blocks
        if line.strip().startswith('```'):
            if not in_code_block:
                in_code_block = True
                code_lines = []
            else:
                in_code_block = False
                # Print code block
                pdf.set_font('DejaVu', '', 8)
                pdf.set_fill_color(240, 240, 240)
                for cl in code_lines:
                    pdf.cell(0, 4, cl, 0, 1, 'L', fill=True)
                pdf.ln(2)
                pdf.set_font('DejaVu', '', 10)
            continue
            
        if in_code_block:
            code_lines.append(line)
            continue
            
        # Handle headers
        if line.startswith('# '):
            pdf.set_font('DejaVu', 'B', 16)
            pdf.cell(0, 10, line[2:], 0, 1, 'L')
            pdf.ln(3)
            pdf.set_font('DejaVu', '', 10)
        elif line.startswith('## '):
            pdf.set_font('DejaVu', 'B', 14)
            pdf.cell(0, 8, line[3:], 0, 1, 'L')
            pdf.ln(2)
            pdf.set_font('DejaVu', '', 10)
        elif line.startswith('### '):
            pdf.set_font('DejaVu', 'B', 12)
            pdf.cell(0, 7, line[4:], 0, 1, 'L')
            pdf.ln(1)
            pdf.set_font('DejaVu', '', 10)
        elif line.startswith('#### '):
            pdf.set_font('DejaVu', 'B', 11)
            pdf.cell(0, 6, line[5:], 0, 1, 'L')
            pdf.ln(1)
            pdf.set_font('DejaVu', '', 10)
        # Handle bold/italic
        elif line.strip().startswith('- ') or line.strip().startswith('* '):
            pdf.set_font('DejaVu', '', 10)
            pdf.cell(10, 5, '', 0, 0)
            pdf.cell(0, 5, line.strip()[2:], 0, 1)
        elif line.strip().startswith('1. ') or line.strip().startswith('2. ') or line.strip().startswith('3. '):
            pdf.set_font('DejaVu', '', 10)
            pdf.cell(10, 5, '', 0, 0)
            pdf.cell(0, 5, line.strip(), 0, 1)
        # Handle tables (simplified)
        elif '|' in line and line.count('|') > 1:
            pdf.set_font('DejaVu', '', 8)
            cells = [c.strip() for c in line.split('|') if c.strip()]
            if cells:
                col_width = 180 / len(cells)
                for cell in cells:
                    pdf.cell(col_width, 5, cell[:30], 1, 0, 'C')
                pdf.ln()
        # Regular text
        elif line.strip():
            pdf.set_font('DejaVu', '', 10)
            # Handle inline formatting
            text = line
            # Simple word wrap
            pdf.multi_cell(0, 5, text)
        else:
            pdf.ln(2)
    
    pdf.output(output_path)
    print(f"Generated: {output_path}")

def main():
    base_dir = Path(__file__).parent.parent
    
    # Files to convert
    files = {
        'AMR_BL_STABILITY_REPORT_TR.md': 'AMR_BL_STABILITY_REPORT_TR.pdf',
        'AMR_BL_STABILITY_REPORT_EN.md': 'AMR_BL_STABILITY_REPORT_EN.pdf',
        'AMR_BL_STABILITY_REPORT.md': 'AMR_BL_STABILITY_REPORT_BILINGUAL.pdf',
    }
    
    for md_file, pdf_file in files.items():
        md_path = base_dir / md_file
        pdf_path = base_dir / pdf_file
        if md_path.exists():
            content = read_markdown(md_path)
            markdown_to_pdf(content, pdf_path)
        else:
            print(f"Warning: {md_file} not found")
    
    print("\nPDF generation complete!")

if __name__ == '__main__':
    main()