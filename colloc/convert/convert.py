import fitz
import os

docs_dir = os.path.join(os.path.dirname(__file__), '..', 'docs')
output_dir = os.path.join(os.path.dirname(__file__), 'converted')
os.makedirs(output_dir, exist_ok=True)

for fname in sorted(os.listdir(docs_dir)):
    if not fname.endswith('.pdf'):
        continue
    fpath = os.path.join(docs_dir, fname)
    doc = fitz.open(fpath)
    print(f"Processing: {fname} ({doc.page_count} pages)")
    full_text = []
    for i, page in enumerate(doc):
        text = page.get_text()
        full_text.append(f"--- Page {i+1} ---\n{text}")
    doc.close()
    out_text = "\n".join(full_text)
    out_name = fname.replace('.pdf', '.txt')
    out_path = os.path.join(output_dir, out_name)
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write(out_text)
    print(f"  -> {out_name} ({len(out_text)} chars)")

print("\nDone. All texts saved in:", output_dir)
