"""
Convert manuscript_draft_v1.md to publication-quality PDF for bioRxiv.
Uses fpdf2 with full Unicode support.
"""
import re
from pathlib import Path
from fpdf import FPDF

PROJECT_ROOT = Path('e:/miniconda3/envs/llama-env/genesis_project')
MD_FILE  = PROJECT_ROOT / 'science_engine/reports/manuscript_draft_v1.md'
OUT_PDF  = PROJECT_ROOT / 'science_engine/reports/manuscript_draft_v1.pdf'


# ── Colour palette ────────────────────────────────────────────────────────────
BLACK  = (0,   0,   0)
DARK   = (30,  30,  30)
GREY   = (90,  90,  90)
LGREY  = (180, 180, 180)
HEAD1  = (20,  60,  120)
HEAD2  = (40,  80,  150)
HEAD3  = (60, 100, 160)
ACCENT = (200, 230, 255)


FONT_DIR = Path('C:/Windows/Fonts')

class ManuscriptPDF(FPDF):
    def __init__(self):
        super().__init__(orientation='P', unit='mm', format='A4')
        self.set_margins(22, 22, 22)
        self.set_auto_page_break(True, margin=22)
        # Register Arial (Unicode-capable TTF)
        self.add_font('Arial',  '',   str(FONT_DIR / 'arial.ttf'))
        self.add_font('Arial',  'B',  str(FONT_DIR / 'arialbd.ttf'))
        self.add_font('Arial',  'I',  str(FONT_DIR / 'ariali.ttf'))
        self.add_font('Arial',  'BI', str(FONT_DIR / 'arialbi.ttf'))
        self.FONT = 'Arial'

    def header(self):
        if self.page_no() > 1:
            self.set_font(self.FONT, 'I', 8)
            self.set_text_color(*GREY)
            self.cell(0, 5, 'ISS Bacillales Defence System Streamlining  \u00b7  ZJY  \u00b7  bioRxiv preprint', align='L')
            self.ln(3)
            self.set_draw_color(*LGREY)
            self.line(self.l_margin, self.get_y(), self.w - self.r_margin, self.get_y())
            self.ln(3)

    def footer(self):
        self.set_y(-15)
        self.set_font(self.FONT, 'I', 8)
        self.set_text_color(*GREY)
        self.cell(0, 10, f'Page {self.page_no()}', align='C')


def clean(text):
    """Remove markdown inline formatting and normalize Unicode for Arial."""
    text = re.sub(r'\*\*(.+?)\*\*', r'\1', text)
    text = re.sub(r'\*(.+?)\*',   r'\1', text)
    text = re.sub(r'_(.+?)_',     r'\1', text)
    text = re.sub(r'`(.+?)`',     r'\1', text)
    # Replace superscript citation markers with bracketed numbers
    sup_map = {'⁰':'0','¹':'1','²':'2','³':'3','⁴':'4',
               '⁵':'5','⁶':'6','⁷':'7','⁸':'8','⁹':'9','⁻':'-',
               '₀':'0','₁':'1','₂':'2','₃':'3','₄':'4',
               '₅':'5','₆':'6','₇':'7','₈':'8','₉':'9'}
    for k, v in sup_map.items():
        text = text.replace(k, v)
    return text.strip()


def inline_bold_chunks(text):
    """Split text into [(is_bold, chunk), ...] for mixed rendering."""
    pattern = r'\*\*(.+?)\*\*'
    parts = []
    last = 0
    for m in re.finditer(pattern, text):
        if m.start() > last:
            parts.append((False, text[last:m.start()]))
        parts.append((True, m.group(1)))
        last = m.end()
    if last < len(text):
        parts.append((False, text[last:]))
    return parts if parts else [(False, text)]


def italic_clean(text):
    """Return (is_italic, cleaned_text) pairs for *italic* inline."""
    pattern = r'\*(.+?)\*'
    parts = []
    last = 0
    for m in re.finditer(pattern, text):
        if m.start() > last:
            parts.append(('', text[last:m.start()]))
        parts.append(('I', m.group(1)))
        last = m.end()
    if last < len(text):
        parts.append(('', text[last:]))
    return parts if parts else [('', text)]


def write_mixed(pdf, text, base_size=10, base_style=''):
    """Write text with inline bold (**)."""
    # First pass: bold chunks
    chunks = inline_bold_chunks(text)
    for is_bold, chunk in chunks:
        # Second pass: italic within each chunk
        style = 'B' if is_bold else base_style
        sub = italic_clean(chunk)
        for it, s in sub:
            if not s:
                continue
            combined = ''
            if 'B' in style and it == 'I':
                combined = 'BI'
            elif 'B' in style:
                combined = 'B'
            elif it == 'I':
                combined = 'I'
            pdf.set_font(pdf.FONT, combined, base_size)
            pdf.write(5, s)
    pdf.set_font(pdf.FONT, base_style, base_size)


def render_table(pdf, rows, col_widths=None):
    """Render a markdown table."""
    if len(rows) < 2:
        return
    header = rows[0]
    body   = rows[2:]  # skip separator row

    usable = pdf.w - pdf.l_margin - pdf.r_margin
    n_cols = len(header)
    if col_widths is None:
        col_widths = [usable / n_cols] * n_cols

    pdf.ln(2)
    pdf.set_draw_color(*LGREY)

    # Header row
    pdf.set_fill_color(*ACCENT)
    pdf.set_font(pdf.FONT, 'B', 8)
    pdf.set_text_color(*HEAD1)
    x0 = pdf.get_x()
    y0 = pdf.get_y()
    for i, cell in enumerate(header):
        pdf.set_xy(x0 + sum(col_widths[:i]), y0)
        pdf.cell(col_widths[i], 6, clean(cell)[:40], border=1, align='C', fill=True)
    pdf.ln(6)

    # Body rows
    pdf.set_font(pdf.FONT, '', 7.5)
    pdf.set_text_color(*DARK)
    for r, row in enumerate(body):
        x0 = pdf.get_x()
        y0 = pdf.get_y()
        fill = r % 2 == 0
        pdf.set_fill_color(248, 248, 252) if fill else pdf.set_fill_color(255, 255, 255)
        for i, cell in enumerate(row):
            pdf.set_xy(x0 + sum(col_widths[:i]), y0)
            txt = clean(cell)
            pdf.cell(col_widths[i], 5.5, txt[:50], border=1, align='C', fill=True)
        pdf.ln(5.5)
    pdf.set_text_color(*DARK)
    pdf.ln(2)


def build_pdf(md_path, out_path):
    text = md_path.read_text(encoding='utf-8')
    lines = text.split('\n')

    pdf = ManuscriptPDF()
    pdf.add_page()

    i = 0
    in_abstract = False
    abstract_lines = []
    table_buffer = []
    in_table = False
    skip_next_hr = False

    while i < len(lines):
        line = lines[i]
        raw  = line.rstrip()

        # ── Horizontal rule ───────────────────────────────────────────────────
        if raw == '---':
            if in_abstract and abstract_lines:
                # Render collected abstract
                pdf.set_fill_color(*ACCENT)
                pdf.set_draw_color(*HEAD2)
                pdf.set_line_width(0.3)
                bx = pdf.l_margin
                by = pdf.get_y()
                bw = pdf.w - pdf.l_margin - pdf.r_margin
                # measure height
                pdf.set_font(pdf.FONT, '', 9)
                full_text = ' '.join(abstract_lines)
                pdf.set_xy(bx + 3, by + 3)
                pdf.multi_cell(bw - 6, 4.5, clean(full_text))
                bh = pdf.get_y() - by + 3
                pdf.rect(bx, by, bw, bh, style='DF')
                # re-render text on top
                pdf.set_xy(bx + 3, by + 3)
                pdf.set_font(pdf.FONT, '', 9)
                pdf.multi_cell(bw - 6, 4.5, clean(full_text))
                pdf.ln(4)
                in_abstract = False
                abstract_lines = []
            else:
                pdf.set_draw_color(*LGREY)
                pdf.set_line_width(0.2)
                pdf.line(pdf.l_margin, pdf.get_y(), pdf.w - pdf.r_margin, pdf.get_y())
                pdf.ln(3)
            i += 1
            continue

        # ── Table detection ───────────────────────────────────────────────────
        if raw.startswith('|'):
            table_buffer.append([c.strip() for c in raw.strip('|').split('|')])
            i += 1
            continue
        else:
            if table_buffer:
                render_table(pdf, table_buffer)
                table_buffer = []

        # ── Heading 1 (title) ─────────────────────────────────────────────────
        if raw.startswith('# ') and not raw.startswith('## '):
            pdf.set_font(pdf.FONT, 'B', 15)
            pdf.set_text_color(*HEAD1)
            pdf.multi_cell(0, 7, clean(raw[2:]))
            pdf.ln(4)
            i += 1
            continue

        # ── Heading 2 ─────────────────────────────────────────────────────────
        if raw.startswith('## '):
            if in_abstract:
                in_abstract = False
            pdf.ln(3)
            pdf.set_font(pdf.FONT, 'B', 12)
            pdf.set_text_color(*HEAD1)
            title = clean(raw[3:])
            pdf.cell(0, 7, title, ln=True)
            # underline
            pdf.set_draw_color(*HEAD2)
            pdf.set_line_width(0.4)
            pdf.line(pdf.l_margin, pdf.get_y(), pdf.l_margin + 80, pdf.get_y())
            pdf.ln(3)
            if title == 'Abstract':
                in_abstract = True
            i += 1
            continue

        # ── Heading 3 ─────────────────────────────────────────────────────────
        if raw.startswith('### '):
            pdf.ln(2)
            pdf.set_font(pdf.FONT, 'BI', 10.5)
            pdf.set_text_color(*HEAD2)
            pdf.cell(0, 6, clean(raw[4:]), ln=True)
            pdf.set_text_color(*DARK)
            pdf.ln(1)
            i += 1
            continue

        # ── Heading 4 ─────────────────────────────────────────────────────────
        if raw.startswith('#### '):
            pdf.set_font(pdf.FONT, 'B', 10)
            pdf.set_text_color(*HEAD3)
            pdf.cell(0, 5, clean(raw[5:]), ln=True)
            pdf.set_text_color(*DARK)
            i += 1
            continue

        # ── Bold metadata lines (Author, Keywords, Running title) ─────────────
        if raw.startswith('**') and raw.endswith('**') and len(raw) > 4:
            pdf.set_font(pdf.FONT, 'B', 9)
            pdf.set_text_color(*GREY)
            pdf.cell(0, 5, clean(raw), ln=True)
            pdf.set_text_color(*DARK)
            i += 1
            continue

        # ── Numbered reference line ───────────────────────────────────────────
        if re.match(r'^\d+\.', raw):
            if in_abstract:
                abstract_lines.append(raw)
            else:
                pdf.set_font(pdf.FONT, '', 8.5)
                pdf.set_text_color(*DARK)
                # indent reference block
                pdf.set_x(pdf.l_margin + 5)
                pdf.multi_cell(pdf.w - pdf.l_margin - pdf.r_margin - 5, 4.5,
                               clean(raw))
            i += 1
            continue

        # ── Bullet list ───────────────────────────────────────────────────────
        if raw.startswith('- ') or raw.startswith('* '):
            content = raw[2:]
            pdf.set_font(pdf.FONT, '', 9.5)
            pdf.set_text_color(*DARK)
            pdf.set_x(pdf.l_margin + 4)
            pdf.cell(4, 5, '\u2022')  # bullet •
            pdf.set_x(pdf.l_margin + 8)
            pdf.multi_cell(pdf.w - pdf.l_margin - pdf.r_margin - 8, 5, clean(content))
            i += 1
            continue

        # ── Empty line ────────────────────────────────────────────────────────
        if raw.strip() == '':
            if in_abstract:
                pass
            else:
                pdf.ln(2)
            i += 1
            continue

        # ── Abstract accumulator ──────────────────────────────────────────────
        if in_abstract:
            abstract_lines.append(raw)
            i += 1
            continue

        # ── Regular paragraph ─────────────────────────────────────────────────
        pdf.set_font(pdf.FONT, '', 10)
        pdf.set_text_color(*DARK)
        # Collect continuation lines
        para = raw
        while i + 1 < len(lines):
            nxt = lines[i + 1].rstrip()
            if (nxt == '' or nxt.startswith('#') or nxt.startswith('|')
                    or nxt == '---' or nxt.startswith('- ')
                    or nxt.startswith('**')):
                break
            para += ' ' + nxt
            i += 1
        pdf.set_x(pdf.l_margin)
        pdf.multi_cell(0, 5, clean(para))
        pdf.ln(1)
        i += 1

    # Flush any remaining table
    if table_buffer:
        render_table(pdf, table_buffer)

    pdf.output(str(out_path))
    print(f'PDF saved: {out_path}')
    print(f'Size: {round(out_path.stat().st_size / 1024)} KB  |  Pages: {pdf.page_no()}')


if __name__ == '__main__':
    build_pdf(MD_FILE, OUT_PDF)
