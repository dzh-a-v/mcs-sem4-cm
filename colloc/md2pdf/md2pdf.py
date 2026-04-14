import os
import sys
import re
import markdown
from xhtml2pdf import pisa


def md_to_html(md_path):
    """Convert a Markdown file to an HTML string."""
    with open(md_path, "r", encoding="utf-8") as f:
        content = f.read()

    md_extensions = {
        "tables": {},
        "fenced_code": {},
        "codehilite": {"guess_lang": False},
        "toc": {"permalink": False},
        "sane_lists": {},
        "attr_list": {},
        "md_in_html": {},
    }

    html_body = markdown.markdown(content, extensions=md_extensions)

    html_doc = f"""<!DOCTYPE html>
<html lang="ru">
<head>
    <meta charset="UTF-8">
    <style>
        @page {{
            size: A4;
            margin: 2cm;
            @frame footer {{
                -pdf-frame-content: footerContent;
                bottom: 2cm;
                margin-left: 2cm;
                margin-right: 2cm;
                height: 1cm;
            }}
        }}
        body {{
            font-family: "Arial", sans-serif;
            font-size: 11pt;
            line-height: 1.5;
            color: #000;
        }}
        h1 {{
            font-size: 16pt;
            text-align: center;
            margin-top: 18pt;
            margin-bottom: 12pt;
            page-break-before: always;
            color: #000;
        }}
        h1:first-of-type {{
            page-break-before: avoid;
        }}
        h2 {{
            font-size: 14pt;
            margin-top: 14pt;
            margin-bottom: 8pt;
            color: #1a1a2e;
            page-break-after: avoid;
        }}
        h3 {{
            font-size: 12pt;
            margin-top: 10pt;
            margin-bottom: 6pt;
            color: #16213e;
            page-break-after: avoid;
        }}
        h4 {{
            font-size: 11pt;
            margin-top: 8pt;
            margin-bottom: 4pt;
            font-style: italic;
            page-break-after: avoid;
        }}
        p {{
            margin: 5pt 0;
            text-align: justify;
        }}
        ul, ol {{
            margin: 5pt 0;
            padding-left: 20pt;
        }}
        li {{
            margin-bottom: 2pt;
        }}
        hr {{
            border: none;
            border-top: 1px solid #999;
            margin: 12pt 0;
        }}
        table {{
            border-collapse: collapse;
            width: 100%;
            margin: 8pt 0;
            font-size: 10pt;
        }}
        th, td {{
            border: 1px solid #333;
            padding: 3pt 5pt;
            text-align: left;
        }}
        th {{
            background-color: #e8e8e8;
            font-weight: bold;
        }}
        code {{
            font-family: "Courier New", monospace;
            font-size: 10pt;
            background-color: #f0f0f0;
            padding: 1pt 3pt;
        }}
        pre {{
            background-color: #f0f0f0;
            padding: 8pt;
            overflow-wrap: break-word;
            white-space: pre-wrap;
            font-size: 9pt;
            margin: 6pt 0;
        }}
        pre code {{
            background-color: transparent;
            padding: 0;
        }}
        strong {{
            font-weight: bold;
        }}
        em {{
            font-style: italic;
        }}
        blockquote {{
            margin-left: 15pt;
            padding-left: 10pt;
            border-left: 3px solid #ccc;
            color: #555;
        }}
    </style>
</head>
<body>
{html_body}
</body>
</html>"""

    return html_doc


def md_to_pdf(md_path, pdf_path):
    """Convert a Markdown file to PDF."""
    html_doc = md_to_html(md_path)

    with open(pdf_path, "w+b") as out_file:
        pisa.CreatePDF(html_doc, dest=out_file, encoding="utf-8")

    print(f"PDF saved to: {pdf_path}")


if __name__ == "__main__":
    base_dir = os.path.dirname(os.path.abspath(__file__))
    project_dir = os.path.dirname(base_dir)

    md_file = os.path.join(project_dir, "notes", "colloc.md")
    pdf_file = os.path.join(project_dir, "notes", "colloc.pdf")

    if len(sys.argv) > 1:
        md_file = sys.argv[1]
    if len(sys.argv) > 2:
        pdf_file = sys.argv[2]

    if not os.path.isfile(md_file):
        print(f"Error: Markdown file not found: {md_file}")
        sys.exit(1)

    md_to_pdf(md_file, pdf_file)
