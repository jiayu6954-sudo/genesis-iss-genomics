"""
Convert manuscript_draft_v1.md → publication PDF for bioRxiv submission.
Features: line numbers, embedded figures, table page-break safety,
          species italics, correct scientific notation, author affiliation.
"""
import re
from pathlib import Path
from fpdf import FPDF
from fpdf.enums import XPos, YPos

PROJECT_ROOT = Path('e:/miniconda3/envs/llama-env/genesis_project')
MD_FILE  = PROJECT_ROOT / 'science_engine/reports/manuscript_draft_v1.md'
FIG_DIR  = PROJECT_ROOT / 'science_engine/figures'
OUT_PDF  = PROJECT_ROOT / 'science_engine/reports/manuscript_draft_v1.pdf'
FONT_DIR = Path('C:/Windows/Fonts')

# ── Colour palette ─────────────────────────────────────────────────────────────
DARK   = (30,  30,  30)
GREY   = (100, 100, 100)
LGREY  = (180, 180, 180)
HEAD1  = (20,  60,  120)
HEAD2  = (40,  80,  150)
HEAD3  = (60, 100, 160)
BLUE_BG= (220, 234, 255)
LN_COL = (160, 160, 160)

# Layout constants
L_MARGIN  = 25    # extra space on left for line numbers
R_MARGIN  = 18
T_MARGIN  = 20
LN_X      = 4     # x-position for line number text
LN_W      = 13    # width of line-number column
LINE_H    = 5.0   # mm per body line
LINE_H_SM = 4.5   # mm for small text


# ── Unicode normalisation ──────────────────────────────────────────────────────
def normalise(text: str) -> str:
    """Convert Unicode special chars to Arial-safe equivalents."""
    # Superscripts → regular digits/signs
    sup = {'⁰':'0','¹':'1','²':'2','³':'3','⁴':'4',
           '⁵':'5','⁶':'6','⁷':'7','⁸':'8','⁹':'9','⁻':'-',
           '₀':'0','₁':'1','₂':'2','₃':'3','₄':'4',
           '₅':'5','₆':'6','₇':'7','₈':'8','₉':'9'}
    for k, v in sup.items():
        text = text.replace(k, v)
    # Dashes / quotes
    text = text.replace('\u2013', '\u2013')  # en-dash — Arial OK
    text = text.replace('\u2014', '\u2014')  # em-dash — Arial OK
    text = text.replace('\u2212', '-')       # minus sign → hyphen
    text = text.replace('\u2032', "'")       # prime
    text = text.replace('\u2019', "'")       # right single quote
    text = text.replace('\u201c', '"')
    text = text.replace('\u201d', '"')
    # Greek letters — keep δ, α, β etc. (Arial supports them)
    # Scientific notation: "10-5" → render cleanly, no extra spaces
    text = re.sub(r'(\d)\s*-\s*(\d)', r'\1-\2', text)  # collapse "10 - 5" → "10-5"
    return text


def clean(text: str) -> str:
    """Strip markdown formatting and normalise."""
    text = re.sub(r'\*\*(.+?)\*\*', r'\1', text)
    text = re.sub(r'\*([^*]+?)\*',  r'\1', text)
    text = re.sub(r'_([^_]+?)_',    r'\1', text)
    text = re.sub(r'`([^`]+?)`',    r'\1', text)
    return normalise(text.strip())


def split_inline(text: str):
    """
    Parse inline markdown into segments: (style, content)
    style: '' | 'B' | 'I' | 'BI'
    Handles **bold**, *italic*, ***bold-italic***.
    """
    segments = []
    pattern = re.compile(r'(\*\*\*(.+?)\*\*\*|\*\*(.+?)\*\*|\*([^*\n]+?)\*)')
    last = 0
    for m in pattern.finditer(text):
        if m.start() > last:
            segments.append(('', normalise(text[last:m.start()])))
        if m.group(2):   # ***bold-italic***
            segments.append(('BI', normalise(m.group(2))))
        elif m.group(3): # **bold**
            segments.append(('B',  normalise(m.group(3))))
        elif m.group(4): # *italic*
            segments.append(('I',  normalise(m.group(4))))
        last = m.end()
    if last < len(text):
        segments.append(('', normalise(text[last:])))
    return segments


# ── PDF class ──────────────────────────────────────────────────────────────────
class ManuscriptPDF(FPDF):

    def __init__(self):
        super().__init__(orientation='P', unit='mm', format='A4')
        self.set_margins(L_MARGIN, T_MARGIN, R_MARGIN)
        self.set_auto_page_break(True, margin=22)
        self.add_font('Arial',  '',   str(FONT_DIR / 'arial.ttf'))
        self.add_font('Arial',  'B',  str(FONT_DIR / 'arialbd.ttf'))
        self.add_font('Arial',  'I',  str(FONT_DIR / 'ariali.ttf'))
        self.add_font('Arial',  'BI', str(FONT_DIR / 'arialbi.ttf'))
        self.F = 'Arial'
        self._line_no   = 0     # running line counter
        self._last_y    = T_MARGIN
        self._ln_active = False  # only number main text pages

    # ── Header / Footer ─────────────────────────────────────────────────────
    def header(self):
        if self.page_no() > 1:
            self.set_font(self.F, 'I', 8)
            self.set_text_color(*GREY)
            self.cell(0, 5,
                'Defence system streamlining in ISS Bacillales  \u00b7  ZJY  \u00b7  bioRxiv',
                align='L')
            self.ln(3)
            self.set_draw_color(*LGREY)
            self.line(self.l_margin, self.get_y(),
                      self.w - self.r_margin, self.get_y())
            self.ln(3)
            self._last_y = self.get_y()

    def footer(self):
        self.set_y(-15)
        self.set_font(self.F, 'I', 8)
        self.set_text_color(*GREY)
        self.cell(0, 10, f'{self.page_no()}', align='C')

    # ── Line numbering ──────────────────────────────────────────────────────
    def _emit_line_numbers(self, y_start, y_end, lh=LINE_H):
        """Draw line numbers in left gutter for the y-span just written."""
        if not self._ln_active or y_end <= y_start:
            return
        n = max(1, round((y_end - y_start) / lh))
        self.set_font(self.F, '', 7)
        self.set_text_color(*LN_COL)
        for j in range(n):
            self._line_no += 1
            if self._line_no % 5 == 0:
                y = y_start + j * lh
                self.set_xy(LN_X, y)
                self.cell(LN_W, lh, str(self._line_no), align='R')
        self.set_text_color(*DARK)

    # ── Text helpers ────────────────────────────────────────────────────────
    def write_inline(self, segments, lh=LINE_H, size=10):
        """Write inline-styled segments on the current line."""
        for style, txt in segments:
            if not txt:
                continue
            self.set_font(self.F, style, size)
            self.write(lh, txt)
        self.set_font(self.F, '', size)

    def para(self, text: str, size=10, lh=LINE_H, indent=0):
        """Render a full paragraph with inline styling + line numbers."""
        if not text.strip():
            return
        segs = split_inline(text)
        self.set_x(self.l_margin + indent)
        y0 = self.get_y()
        self.write_inline(segs, lh, size)
        self.ln(lh)
        self._emit_line_numbers(y0, self.get_y(), lh)

    def heading(self, text: str, level=2):
        """Render a section heading."""
        if level == 1:
            self.ln(2)
            self.set_font(self.F, 'B', 15)
            self.set_text_color(*HEAD1)
            y0 = self.get_y()
            self.multi_cell(0, 7, clean(text))
            self.ln(3)
        elif level == 2:
            self.ln(4)
            self.set_font(self.F, 'B', 12)
            self.set_text_color(*HEAD1)
            y0 = self.get_y()
            self.cell(0, 7, clean(text),
                      new_x=XPos.LMARGIN, new_y=YPos.NEXT)
            self.set_draw_color(*HEAD2)
            self.set_line_width(0.4)
            self.line(self.l_margin, self.get_y(),
                      self.l_margin + 90, self.get_y())
            self.ln(3)
        elif level == 3:
            self.ln(2)
            self.set_font(self.F, 'BI', 10.5)
            self.set_text_color(*HEAD2)
            y0 = self.get_y()
            self.cell(0, 6, clean(text),
                      new_x=XPos.LMARGIN, new_y=YPos.NEXT)
            self.ln(1)
        elif level == 4:
            self.set_font(self.F, 'B', 10)
            self.set_text_color(*HEAD3)
            self.cell(0, 5, clean(text),
                      new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.set_text_color(*DARK)

    # ── Table ───────────────────────────────────────────────────────────────
    def table(self, rows):
        """Render a markdown table. Force page break if it won't fit."""
        if len(rows) < 3:
            return
        header  = rows[0]
        body    = [r for r in rows[2:] if r and any(c.strip() for c in r)]
        n_cols  = len(header)
        if n_cols == 0:
            return

        usable  = self.w - self.l_margin - self.r_margin
        col_w   = [usable / n_cols] * n_cols

        # Wider first column for most tables
        if n_cols >= 4:
            col_w[0] = usable * 0.28
            rest = (usable - col_w[0]) / (n_cols - 1)
            for k in range(1, n_cols):
                col_w[k] = rest

        row_h  = 5.5
        th     = 6 + (len(body)) * row_h + 6  # estimated table height
        if self.get_y() + th > self.h - 25:
            self.add_page()

        self.ln(2)
        self.set_draw_color(*LGREY)
        self.set_line_width(0.2)

        # Header
        self.set_fill_color(*BLUE_BG)
        self.set_font(self.F, 'B', 8)
        self.set_text_color(*HEAD1)
        x0, y0 = self.get_x(), self.get_y()
        for k, cell in enumerate(header):
            self.set_xy(x0 + sum(col_w[:k]), y0)
            self.cell(col_w[k], 6, clean(cell)[:45],
                      border=1, align='C', fill=True)
        self.set_xy(x0, y0 + 6)

        # Body
        self.set_font(self.F, '', 7.5)
        self.set_text_color(*DARK)
        for ridx, row in enumerate(body):
            # pad / trim row to n_cols
            row = (row + [''] * n_cols)[:n_cols]
            fill_c = (248, 248, 252) if ridx % 2 == 0 else (255, 255, 255)
            self.set_fill_color(*fill_c)
            x0, y0 = self.get_x(), self.get_y()
            for k, cell in enumerate(row):
                self.set_xy(x0 + sum(col_w[:k]), y0)
                txt = clean(cell)[:55]
                self.cell(col_w[k], row_h, txt,
                          border=1, align='C', fill=True)
            self.set_xy(x0, y0 + row_h)
        self.ln(4)

    # ── Figure embed ────────────────────────────────────────────────────────
    def embed_figure(self, png_path: Path, caption: str, label: str):
        """Embed a PNG figure with caption on its own page."""
        self.add_page()
        self._ln_active = False
        # Label
        self.set_font(self.F, 'B', 10)
        self.set_text_color(*HEAD1)
        self.cell(0, 6, label, new_x=XPos.LMARGIN, new_y=YPos.NEXT)
        self.ln(2)
        # Image
        if png_path.exists():
            avail_w = self.w - self.l_margin - self.r_margin
            avail_h = self.h - self.get_y() - 35
            self.image(str(png_path), x=self.l_margin, w=avail_w,
                       h=min(avail_h, 120))
            self.ln(4)
        else:
            self.set_font(self.F, 'I', 9)
            self.set_text_color(*GREY)
            self.cell(0, 6, f'[Figure file not found: {png_path.name}]',
                      new_x=XPos.LMARGIN, new_y=YPos.NEXT)
            self.ln(2)
        # Caption
        self.set_font(self.F, '', 9)
        self.set_text_color(*DARK)
        segs = split_inline(caption)
        self.write_inline(segs, LINE_H_SM, 9)
        self.ln(LINE_H_SM)
        self._ln_active = True


# ── Main builder ───────────────────────────────────────────────────────────────
def build_pdf(md_path: Path, out_path: Path):
    text  = md_path.read_text(encoding='utf-8')
    lines = text.split('\n')

    pdf = ManuscriptPDF()
    pdf.add_page()
    pdf._ln_active = True

    i          = 0
    in_abstract= False
    abs_buf    = []
    tbl_buf    = []
    author_done= False
    fig_done   = False   # have we inserted figure pages yet?

    # Figure captions from the manuscript
    fig_data = [
        (FIG_DIR / 'Figure1.png',
         'Figure 1. Species-stratified analysis of CRISPR-Cas and IS element '
         'distributions in ISS versus ground Bacillales. '
         '(a) Proportion of strains with at least one Cas-encoding CDS, '
         'stratified by species and environment. Only *P. polymyxa* shows a '
         'significant ISS-vs-ground difference (Fisher p = 0.0095; BH q = 0.025). '
         '(b) Boxplots of IS element density (IS/1,000 CDS) by species and '
         'environment. ISS = orange, Ground = blue.',
         'Figure 1'),
        (FIG_DIR / 'Figure2.png',
         'Figure 2. Co-loss of CRISPR-Cas and IS elements in *Paenibacillus '
         'polymyxa* and metagenomics validation. '
         '(a) Scatter plot of IS element density vs CRISPR-Cas status in all '
         '*P. polymyxa* strains (p = 0.003, Cliff delta = -0.75). '
         '(b) CRISPR-associated KO density in GLDS-224 ISS debris vs ground '
         'debris (ratio 0.503-fold).',
         'Figure 2'),
        (FIG_DIR / 'FigureS1.png',
         'Supplementary Figure S1. Assembly quality assessment. '
         '(a) N50 vs IS density scatter across all 85 genomes (r = +0.428). '
         '(b) IS/1,000 CDS vs IS/Mbp normalisation concordance by species.',
         'Supplementary Figure S1'),
    ]

    while i < len(lines):
        raw = lines[i].rstrip()

        # ── Flush table buffer ─────────────────────────────────────────────
        if not raw.startswith('|') and tbl_buf:
            pdf.table(tbl_buf)
            tbl_buf = []

        # ── Horizontal rule ────────────────────────────────────────────────
        if raw == '---':
            if in_abstract and abs_buf:
                # Render abstract box
                full = normalise(' '.join(abs_buf))
                segs = split_inline(' '.join(abs_buf))
                bx, by = pdf.l_margin, pdf.get_y()
                bw = pdf.w - pdf.l_margin - pdf.r_margin
                # Estimate height from character count (no ghost render)
                pdf.set_font(pdf.F, '', 9)
                avg_cpl = int((bw - 6) / 1.85)  # ~chars per line at 9pt Arial
                n_lines = max(1, -(-len(full) // avg_cpl))  # ceil division
                bh = n_lines * LINE_H_SM + 8
                # Draw background rect first
                pdf.set_fill_color(*BLUE_BG)
                pdf.set_draw_color(*HEAD2)
                pdf.set_line_width(0.3)
                pdf.rect(bx, by, bw, bh, style='DF')
                # Render text once on top
                pdf.set_xy(bx + 3, by + 3)
                pdf.set_font(pdf.F, '', 9)
                pdf.multi_cell(bw - 6, LINE_H_SM, full)
                pdf.ln(4)
                in_abstract = False
                abs_buf = []
            else:
                pdf.set_draw_color(*LGREY)
                pdf.set_line_width(0.2)
                pdf.line(pdf.l_margin, pdf.get_y(),
                         pdf.w - pdf.r_margin, pdf.get_y())
                pdf.ln(3)
            i += 1
            continue

        # ── Insert figures before Figure Legends section ───────────────────
        if not fig_done and raw.startswith('## Figure Legends'):
            for png, cap, lbl in fig_data:
                pdf.embed_figure(png, cap, lbl)
            fig_done = True
            # Fall through to also render the heading

        # ── Table rows ─────────────────────────────────────────────────────
        if raw.startswith('|'):
            tbl_buf.append([c.strip() for c in raw.strip('|').split('|')])
            i += 1
            continue

        # ── Headings ───────────────────────────────────────────────────────
        if raw.startswith('#### '):
            pdf.heading(raw[5:], 4)
            i += 1
            continue
        if raw.startswith('### '):
            pdf.heading(raw[4:], 3)
            i += 1
            continue
        if raw.startswith('## '):
            title = raw[3:]
            pdf.heading(title, 2)
            if title.strip() == 'Abstract':
                in_abstract = True
            i += 1
            continue
        if raw.startswith('# ') and not raw.startswith('## '):
            pdf.heading(raw[2:], 1)
            i += 1
            continue

        # ── Author / Correspondence metadata lines ─────────────────────────
        if raw.startswith('**Author:**'):
            if not author_done:
                pdf.set_font(pdf.F, 'B', 10)
                pdf.set_text_color(*DARK)
                pdf.cell(0, 5, 'Author: ZJY',
                         new_x=XPos.LMARGIN, new_y=YPos.NEXT)
                pdf.set_font(pdf.F, 'I', 9)
                pdf.set_text_color(*GREY)
                pdf.cell(0, 5, 'Affiliation: China',
                         new_x=XPos.LMARGIN, new_y=YPos.NEXT)
                pdf.set_text_color(*DARK)
                author_done = True
            i += 1
            continue
        if raw.startswith('**Correspondence:**'):
            pdf.set_font(pdf.F, '', 9)
            pdf.set_text_color(*GREY)
            pdf.cell(0, 5, 'Correspondence: jiayu6954@gmail.com',
                     new_x=XPos.LMARGIN, new_y=YPos.NEXT)
            pdf.set_text_color(*DARK)
            i += 1
            continue

        # ── Bold-only lines (Running title, Keywords) ──────────────────────
        m_kw = re.match(r'^\*\*([^*]+?)\*\*\s*(.*)', raw)
        if m_kw and raw.endswith('**') is False:
            # e.g. **Running title:** Defence system streamlining...
            if any(kw in raw for kw in ['Running title', 'Keywords']):
                pdf.set_font(pdf.F, 'B', 9)
                pdf.set_text_color(*GREY)
                pdf.cell(0, 5, clean(raw),
                         new_x=XPos.LMARGIN, new_y=YPos.NEXT)
                pdf.set_text_color(*DARK)
                i += 1
                continue

        # ── Numbered reference ─────────────────────────────────────────────
        if re.match(r'^\d+\.', raw):
            if in_abstract:
                abs_buf.append(raw)
            else:
                pdf.set_font(pdf.F, '', 8.5)
                pdf.set_x(pdf.l_margin + 4)
                y0 = pdf.get_y()
                segs = split_inline(raw)
                pdf.write_inline(segs, LINE_H_SM, 8.5)
                pdf.ln(LINE_H_SM)
                pdf._emit_line_numbers(y0, pdf.get_y(), LINE_H_SM)
            i += 1
            continue

        # ── Bullet ─────────────────────────────────────────────────────────
        if raw.startswith('- ') or raw.startswith('* '):
            content = raw[2:]
            pdf.set_font(pdf.F, '', 9.5)
            pdf.set_x(pdf.l_margin + 2)
            pdf.cell(5, LINE_H, '\u2022')
            pdf.set_x(pdf.l_margin + 7)
            y0 = pdf.get_y()
            segs = split_inline(content)
            pdf.write_inline(segs, LINE_H, 9.5)
            pdf.ln(LINE_H)
            pdf._emit_line_numbers(y0, pdf.get_y())
            i += 1
            continue

        # ── Empty line ─────────────────────────────────────────────────────
        if raw.strip() == '':
            if in_abstract:
                pass
            else:
                pdf.ln(2)
            i += 1
            continue

        # ── Abstract accumulator ───────────────────────────────────────────
        if in_abstract:
            abs_buf.append(raw)
            i += 1
            continue

        # ── Regular paragraph (collect continuation lines) ─────────────────
        para = raw
        while i + 1 < len(lines):
            nxt = lines[i + 1].rstrip()
            if (not nxt or nxt.startswith('#') or nxt.startswith('|')
                    or nxt == '---' or nxt.startswith('- ')
                    or nxt.startswith('* ') or re.match(r'^\d+\.', nxt)
                    or nxt.startswith('**')):
                break
            para += ' ' + nxt
            i += 1
        pdf.para(para)
        i += 1

    # Flush remaining table
    if tbl_buf:
        pdf.table(tbl_buf)

    pdf.output(str(out_path))
    size_kb = round(out_path.stat().st_size / 1024)
    print(f'PDF saved  : {out_path}')
    print(f'Size       : {size_kb} KB  |  Pages: {pdf.page_no()}')
    print(f'Line count : {pdf._line_no}')


if __name__ == '__main__':
    build_pdf(MD_FILE, OUT_PDF)
