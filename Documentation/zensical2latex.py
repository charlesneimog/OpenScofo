#!/usr/bin/env python3
"""Convert a Zensical/MkDocs-style project into a small LaTeX project."""

from __future__ import annotations

import argparse
import hashlib
import html
import os
import re
import shutil
import subprocess
import tempfile
import tomllib
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urlparse


LATEX_SPECIALS = {
    "\\": r"\textbackslash{}",
    "&": r"\&",
    "%": r"\%",
    "$": r"\$",
    "#": r"\#",
    "_": r"\_",
    "{": r"\{",
    "}": r"\}",
    "~": r"\textasciitilde{}",
    "^": r"\textasciicircum{}",
    "♯": r"\(\sharp\)",
    "♭": r"\(\flat\)",
    "…": r"\ldots{}",
    "×": r"\(\times\)",
}

HEADING_COMMANDS = {
    1: "chapter",
    2: "section",
    3: "subsection",
    4: "subsubsection",
    5: "paragraph",
    6: "subparagraph",
}

RAW_LATEX_PREFIX = "@@RAW_LATEX_BLOCK@@"
MERMAID_INIT = (
    '%%{init: {"theme": "base", "themeVariables": {"primaryTextColor": "#000000", '
    '"secondaryTextColor": "#000000", "tertiaryTextColor": "#000000", "lineColor": "#000000", '
    '"primaryBorderColor": "#000000", "secondaryBorderColor": "#000000", "tertiaryBorderColor": "#000000"}, '
    '"flowchart": {"htmlLabels": false, "curve": "linear"}}}%%'
)
CURRENT_MARKDOWN_DIR: Path | None = None
CURRENT_MARKDOWN_FILE: Path | None = None
CURRENT_DOCS_DIR: Path | None = None
CURRENT_OUTPUT_DIR: Path | None = None
CURRENT_CHAPTER_LINKS: dict[Path, tuple[str, str]] = {}
CURRENT_FOOTNOTES: dict[str, str] = {}


@dataclass(frozen=True)
class Page:
    source: Path
    title: str
    tex_name: str
    section_path: tuple[str, ...]


@dataclass(frozen=True)
class ExternalLink:
    title: str
    url: str
    section_path: tuple[str, ...]


NavEntry = Page | ExternalLink


def latex_escape(text: str) -> str:
    return "".join(LATEX_SPECIALS.get(char, char) for char in text)


def strip_markdown_formatting(text: str) -> str:
    text = re.sub(r"`([^`]+)`", r"\1", text)
    text = re.sub(r"\*\*([^*]+)\*\*", r"\1", text)
    text = re.sub(r"__([^_]+)__", r"\1", text)
    text = re.sub(r"\*([^*]+)\*", r"\1", text)
    text = re.sub(r"_([^_]+)_", r"\1", text)
    text = re.sub(r":[-+a-zA-Z0-9_]+:", "", text)
    return html.unescape(re.sub(r"<[^>]+>", "", text)).strip()


def slug_for_path(path: Path) -> str:
    stem = "_".join(path.with_suffix("").parts)
    return re.sub(r"[^A-Za-z0-9_]+", "_", stem).strip("_").lower() or "index"


def title_from_markdown(path: Path) -> str:
    try:
        text = path.read_text(encoding="utf-8")
    except FileNotFoundError:
        return path.stem.replace("-", " ").title()

    for line in text.splitlines():
        match = re.match(r"^#\s+(.+?)\s*$", line)
        if match:
            return strip_markdown_formatting(match.group(1))
    return path.stem.replace("-", " ").title()


def is_url(value: str) -> bool:
    parsed = urlparse(value)
    return parsed.scheme in {"http", "https"}


def load_project(toml_path: Path) -> dict:
    with toml_path.open("rb") as handle:
        return tomllib.load(handle)["project"]


def flatten_nav(
    nav: list,
    docs_dir: Path,
    section_path: tuple[str, ...] = (),
) -> list[NavEntry]:
    entries: list[NavEntry] = []
    for item in nav:
        if isinstance(item, str):
            if is_url(item):
                entries.append(ExternalLink(item, item, section_path))
                continue
            source = docs_dir / item
            entries.append(Page(source, title_from_markdown(source), f"{slug_for_path(Path(item))}.tex", section_path))
        elif isinstance(item, dict):
            for title, children in item.items():
                if isinstance(children, str):
                    if is_url(children):
                        entries.append(ExternalLink(title, children, section_path))
                    else:
                        source = docs_dir / children
                        entries.append(
                            Page(source, title_from_markdown(source), f"{slug_for_path(Path(children))}.tex", section_path + (title,))
                        )
                elif isinstance(children, list):
                    entries.extend(flatten_nav(children, docs_dir, section_path + (title,)))
        else:
            raise TypeError(f"Unsupported nav entry: {item!r}")
    return entries


def strip_front_matter(markdown: str) -> str:
    if markdown.startswith("---\n"):
        end = markdown.find("\n---", 4)
        if end != -1:
            return markdown[end + 4 :].lstrip()
    return markdown


def remove_latex_false_html(markdown: str) -> str:
    latex_false_attr = r"\blatex\s*=\s*(?:\"false\"|'false'|false)(?=\s|/|>)"
    tag_name = r"[A-Za-z][A-Za-z0-9:-]*"
    block_pattern = re.compile(
        rf"<(?P<tag>{tag_name})\b(?=[^>]*{latex_false_attr})[^>]*>.*?</(?P=tag)\s*>",
        flags=re.IGNORECASE | re.DOTALL,
    )
    previous = None
    while previous != markdown:
        previous = markdown
        markdown = block_pattern.sub("", markdown)
    return re.sub(rf"<{tag_name}\b(?=[^>]*{latex_false_attr})[^>]*?/?>", "", markdown, flags=re.IGNORECASE)


def path_for_latex(path: Path) -> str:
    if CURRENT_OUTPUT_DIR is None:
        return path.as_posix()
    return Path(os.path.relpath(path, CURRENT_OUTPUT_DIR)).as_posix()


def render_svg(svg_path: Path) -> str | None:
    if CURRENT_OUTPUT_DIR is None or shutil.which("rsvg-convert") is None:
        return None
    generated_dir = CURRENT_OUTPUT_DIR / "generated-images"
    generated_dir.mkdir(parents=True, exist_ok=True)
    try:
        if CURRENT_DOCS_DIR is not None:
            relative = svg_path.resolve().relative_to(CURRENT_DOCS_DIR)
        else:
            relative = Path(svg_path.name)
    except ValueError:
        relative = Path(svg_path.name)
    pdf_path = generated_dir / relative.with_suffix(".pdf")
    pdf_path.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["rsvg-convert", "-f", "pdf", "-o", str(pdf_path), str(svg_path)],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    return path_for_latex(pdf_path)


def render_mermaid(source: str) -> str | None:
    if CURRENT_OUTPUT_DIR is None or shutil.which("mmdc") is None:
        return None
    generated_dir = CURRENT_OUTPUT_DIR / "generated-images" / "mermaid"
    generated_dir.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256(source.encode("utf-8")).hexdigest()[:16]
    mermaid_path = generated_dir / f"{digest}.mmd"
    png_path = generated_dir / f"{digest}.png"
    if png_path.exists():
        return path_for_latex(png_path)

    puppeteer_config = generated_dir / "puppeteer-config.json"
    puppeteer_config.write_text(
        '{"executablePath": "/usr/bin/chromium", "args": ["--no-sandbox", "--disable-setuid-sandbox"]}\n',
        encoding="utf-8",
    )
    mermaid_path.write_text(source, encoding="utf-8")
    result = subprocess.run(
        [
            "mmdc",
            "-i",
            str(mermaid_path),
            "-o",
            str(png_path),
            "-b",
            "white",
            "-s",
            "2",
            "-p",
            str(puppeteer_config),
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    if result.returncode != 0 or not png_path.exists():
        return None
    return path_for_latex(png_path)


def prepare_mermaid_source(source: str) -> str:
    lines = [line for line in source.splitlines() if not re.match(r"^\s*%%\{init:.*\}%%\s*$", line)]
    return "\n".join([MERMAID_INIT, *lines]).strip()


def convert_mermaid(source_lines: list[str]) -> str:
    source = prepare_mermaid_source("\n".join(source_lines).strip())
    rendered = render_mermaid(source)
    if rendered is not None:
        return "\n".join(
            [
                r"\begin{center}",
                rf"\includegraphics[width=0.9\linewidth,height=0.62\textheight,keepaspectratio]{{{latex_escape(rendered)}}}",
                r"\end{center}",
            ]
        )
    return "\n".join([r"\begin{lstlisting}", *source_lines, r"\end{lstlisting}"])


def convert_code_block(language: str, source_lines: list[str]) -> str:
    if language.lower() == "mermaid":
        return convert_mermaid(source_lines)
    if language.lower() == "openscofo":
        options = (
            "listing only, breakable, colback=gray!5, colframe=black!35, "
            "boxrule=0.4pt, arc=1mm, left=1mm, right=1mm, top=1mm, bottom=1mm, "
            r"listing options={basicstyle=\ttfamily\small,breaklines=true,columns=fullflexible,keepspaces=true}"
        )
        return "\n".join(
            [
                rf"\begin{{tcblisting}}{{{options}}}",
                *(sanitize_code_line(line) for line in source_lines),
                r"\end{tcblisting}",
            ]
        )
    return "\n".join([r"\begin{lstlisting}", *(sanitize_code_line(line) for line in source_lines), r"\end{lstlisting}"])


def download_remote_image(url: str) -> str | None:
    suffix = Path(urlparse(url).path).suffix.lower()
    if suffix not in {".png", ".jpg", ".jpeg", ".pdf"}:
        return None
    tmp_dir = Path(tempfile.gettempdir()) / "zensical2latex-images"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256(url.encode("utf-8")).hexdigest()[:16]
    image_path = tmp_dir / f"{digest}{suffix}"
    if not image_path.exists():
        with urllib.request.urlopen(url, timeout=30) as response:
            image_path.write_bytes(response.read())
    return path_for_latex(image_path)


def resolve_local_markdown_link(target: str) -> tuple[str, str] | None:
    if CURRENT_MARKDOWN_DIR is None:
        return None
    target = target.split("{", 1)[0].strip()
    target_path = target.split("#", 1)[0]
    if not target_path or Path(target_path).suffix.lower() != ".md":
        return None
    resolved = (CURRENT_MARKDOWN_DIR / target_path).resolve()
    return CURRENT_CHAPTER_LINKS.get(resolved)


def convert_image(alt: str, target: str) -> str:
    target = target.split("#", 1)[0]
    if not target:
        return ""
    resolved_path: Path | None = None
    if not is_url(target) and CURRENT_MARKDOWN_DIR is not None and CURRENT_DOCS_DIR is not None:
        resolved_path = (CURRENT_MARKDOWN_DIR / target).resolve()
        try:
            target = path_for_latex(resolved_path)
        except ValueError:
            target = target.lstrip("./")
    caption = latex_escape(strip_markdown_formatting(alt))
    if is_url(target):
        downloaded = download_remote_image(target)
        if downloaded is None:
            body = rf"\url{{{target}}}"
        else:
            body = rf"\includegraphics[width=0.8\linewidth]{{{latex_escape(downloaded)}}}"
    elif target.lower().endswith(".svg") and resolved_path is not None:
        rendered = render_svg(resolved_path)
        if rendered is None:
            body = rf"\fbox{{SVG image: \texttt{{{latex_escape(target)}}}}}"
        else:
            body = rf"\includegraphics[width=0.8\linewidth]{{{latex_escape(rendered)}}}"
    else:
        body = rf"\includegraphics[width=0.8\linewidth]{{{latex_escape(target)}}}"
    if caption:
        return "\n".join([r"\begin{figure}[H]", r"\centering", body, rf"\caption{{{caption}}}", r"\end{figure}"])
    return "\n".join([r"\begin{figure}[H]", r"\centering", body, r"\end{figure}"])


def convert_composition_card(image_target: str, link_target: str) -> str:
    figure = convert_image("", image_target).splitlines()
    caption = rf"\caption{{\href{{{link_target}}}{{\texttt{{{latex_escape(link_target)}}}}}}}"
    if figure and figure[-1] == r"\end{figure}":
        figure.insert(-1, caption)
    return "\n".join(figure)


def convert_inline(text: str) -> str:
    placeholders: list[str] = []

    def stash(value: str) -> str:
        placeholders.append(value)
        return f"@@PLACEHOLDER{len(placeholders) - 1}@@"

    def stash_math(match: re.Match[str]) -> str:
        return stash(match.group(0))

    text = html.unescape(text)
    text = re.sub(r"\$\$.*?\$\$", stash_math, text, flags=re.DOTALL)
    text = re.sub(r"\\\[.*?\\\]", stash_math, text, flags=re.DOTALL)
    text = re.sub(r"\\\(.*?\\\)", stash_math, text, flags=re.DOTALL)
    text = re.sub(r"(?<!\\)\$(?!\s)(.+?)(?<!\s)(?<!\\)\$", stash_math, text)

    def code_repl(match: re.Match[str]) -> str:
        return stash(rf"\texttt{{{latex_escape(match.group(1))}}}")

    def bold_repl(match: re.Match[str]) -> str:
        return stash(rf"\textbf{{{convert_inline(match.group(1))}}}")

    def italic_repl(match: re.Match[str]) -> str:
        return stash(rf"\textit{{{convert_inline(match.group(1))}}}")

    def image_repl(match: re.Match[str]) -> str:
        return stash(convert_image(match.group(1), match.group(2)))

    def footnote_repl(match: re.Match[str]) -> str:
        footnote = CURRENT_FOOTNOTES.get(match.group(1))
        if footnote is None:
            return match.group(0)
        return stash(rf"\footnote{{{convert_inline(footnote)}}}")

    def link_repl(match: re.Match[str]) -> str:
        label = convert_inline(match.group(1))
        target = match.group(2).split("{", 1)[0].strip()
        if is_url(target):
            return stash(rf"\href{{{target}}}{{{label}}}")
        local_chapter = resolve_local_markdown_link(target)
        if local_chapter is not None:
            chapter_label, chapter_title = local_chapter
            if not label:
                label = latex_escape(chapter_title)
            return stash(rf"\hyperref[{chapter_label}]{{{label}}}")
        return stash(rf"{label} (\texttt{{{latex_escape(target)}}})")

    text = re.sub(r"!\[([^\]]*)\]\(([^)]+)\)(?:\{[^}]*\})?", image_repl, text)
    text = re.sub(r"\[([^\]]*)\]\(([^)]+)\)(?:\{[^}]*\})?", link_repl, text)
    text = re.sub(r"\[\^([^\]]+)\]", footnote_repl, text)
    text = re.sub(r"`([^`]+)`", code_repl, text)
    text = re.sub(r"\*\*(.+?)\*\*", bold_repl, text)
    text = re.sub(r"__(.+?)__", bold_repl, text)
    text = re.sub(r"(?<!\*)\*([^*\n]+)\*(?!\*)", italic_repl, text)
    text = re.sub(r"(?<!_)_([^_\n]+)_(?!_)", italic_repl, text)
    text = re.sub(r":[-+a-zA-Z0-9_]+:", "", text)

    escaped = latex_escape(text)
    changed = True
    while changed:
        changed = False
        for index, value in enumerate(placeholders):
            escaped_placeholder = latex_escape(f"@@PLACEHOLDER{index}@@")
            raw_placeholder = f"@@PLACEHOLDER{index}@@"
            if escaped_placeholder in escaped or raw_placeholder in escaped:
                escaped = escaped.replace(escaped_placeholder, value).replace(raw_placeholder, value)
                changed = True
    return escaped


def clean_html_line(line: str) -> str:
    line = re.sub(r"<release\b[^>]*>(.*?)</release>", r"\1", line, flags=re.IGNORECASE)
    score = re.search(r"<score\b([^>]*)></score>", line, re.IGNORECASE)
    if score:
        attrs = strip_markdown_formatting(score.group(1)).strip()
        return f"Unsupported score preview: {attrs}"
    line = re.sub(r"<br\s*/?>", r"\\", line, flags=re.IGNORECASE)
    line = re.sub(r"</?hr\s*/?>", "---", line, flags=re.IGNORECASE)
    line = re.sub(r"</?(?:div|p|span|center|style|a|h[1-6])\b[^>]*>", "", line, flags=re.IGNORECASE)
    line = re.sub(r"<i>(.*?)</i>", r"*\1*", line, flags=re.IGNORECASE)
    line = re.sub(r"<em>(.*?)</em>", r"*\1*", line, flags=re.IGNORECASE)
    line = re.sub(r"<b>(.*?)</b>", r"**\1**", line, flags=re.IGNORECASE)
    line = re.sub(r"<strong>(.*?)</strong>", r"**\1**", line, flags=re.IGNORECASE)
    return line


def score_event_image(score_name: str) -> str | None:
    if CURRENT_DOCS_DIR is None or CURRENT_MARKDOWN_FILE is None:
        return None
    try:
        current_doc = CURRENT_MARKDOWN_FILE.resolve().relative_to(CURRENT_DOCS_DIR.resolve())
    except ValueError:
        return None
    if current_doc != Path("score/events.md"):
        return None

    image_name = re.sub(r"[^a-z0-9]+", "_", score_name.lower()).strip("_")
    if not image_name:
        return None
    image_path = CURRENT_DOCS_DIR / "assets" / "events" / f"{image_name}.png"
    if not image_path.exists():
        return None
    return convert_image("", str(image_path))


def sanitize_code_line(line: str) -> str:
    replacements = {
        "├": "|",
        "└": "`",
        "│": "|",
        "─": "-",
        "“": '"',
        "”": '"',
        "‘": "'",
        "’": "'",
        "→": "->",
    }
    return "".join(replacements.get(char, char) for char in line)


def convert_admonition(kind: str, title: str, body: list[str]) -> str:
    heading = title or kind.title()
    converted_body = convert_markdown_lines(body).strip()
    if not converted_body:
        converted_body = r"\vspace{0.1em}"
    return RAW_LATEX_PREFIX + "\n".join(
        [
            rf"\begin{{tcolorbox}}[title={{{convert_inline(heading)}}}, colback=gray!5, colframe=black!35]",
            converted_body,
            r"\end{tcolorbox}",
            "",
        ]
    )


def split_table_row(line: str) -> list[str]:
    line = line.strip()
    if line.startswith("|"):
        line = line[1:]
    if line.endswith("|"):
        line = line[:-1]
    return [cell.strip() for cell in line.split("|")]


def is_table_separator(line: str) -> bool:
    cells = split_table_row(line)
    if not cells:
        return False
    return all(re.match(r"^:?-{3,}:?$", cell.strip()) for cell in cells)


def convert_markdown_table(rows: list[str]) -> str:
    header = split_table_row(rows[0])
    body_rows = [split_table_row(row) for row in rows[2:]]
    column_count = len(header)
    if column_count == 0:
        return ""

    if column_count == 1:
        columns = "|X|"
    else:
        columns = "|" + "|".join(["l"] + ["X"] * (column_count - 1)) + "|"

    def normalize_cells(cells: list[str]) -> list[str]:
        padded = cells[:column_count] + [""] * max(0, column_count - len(cells))
        return padded[:column_count]

    lines = [
        RAW_LATEX_PREFIX + r"\begin{center}",
        rf"\begin{{tabularx}}{{\linewidth}}{{{columns}}}",
        r"\hline",
        " & ".join(rf"\textbf{{{convert_inline(cell)}}}" for cell in normalize_cells(header)) + r" \\",
        r"\hline",
    ]
    for row in body_rows:
        lines.append(" & ".join(convert_inline(cell) for cell in normalize_cells(row)) + r" \\")
        lines.append(r"\hline")
    lines.extend([r"\end{tabularx}", r"\end{center}", ""])
    return "\n".join(lines)


def normalize_tables(lines: list[str]) -> list[str]:
    output: list[str] = []
    index = 0
    in_code = False

    while index < len(lines):
        line = lines[index]
        if re.match(r"^```\s*", line):
            in_code = not in_code
            output.append(line)
            index += 1
            continue

        if not in_code and index + 1 < len(lines) and "|" in line and is_table_separator(lines[index + 1]):
            table_lines = [line, lines[index + 1]]
            index += 2
            while index < len(lines) and "|" in lines[index].strip() and lines[index].strip():
                table_lines.append(lines[index])
                index += 1
            output.append(convert_markdown_table(table_lines))
            continue

        output.append(line)
        index += 1

    return output


def normalize_admonitions(lines: list[str]) -> list[str]:
    output: list[str] = []
    index = 0
    while index < len(lines):
        match = re.match(r"^([!?]{3})\s+([A-Za-z0-9_-]+)(?:\s+\"([^\"]+)\")?\s*$", lines[index])
        if not match:
            output.append(lines[index])
            index += 1
            continue

        kind = match.group(2)
        title = match.group(3) or kind.title()
        index += 1
        body: list[str] = []
        while index < len(lines):
            line = lines[index]
            if line.startswith("    ") or line.startswith("\t") or not line.strip():
                body.append(line[4:] if line.startswith("    ") else line.lstrip("\t"))
                index += 1
                continue
            break
        output.append(convert_admonition(kind, title, body))
    return output


def convert_markdown_lines(lines: list[str]) -> str:
    lines = normalize_admonitions(lines)
    lines = normalize_tables(lines)
    output: list[str] = []
    in_code = False
    code_language = ""
    code_lines: list[str] = []
    in_style = False
    list_stack: list[str] = []
    current_heading = ""

    def close_lists() -> None:
        while list_stack:
            output.append(rf"\end{{{list_stack.pop()}}}")

    for raw_line in lines:
        line = raw_line.rstrip()

        if line.lstrip().startswith(RAW_LATEX_PREFIX):
            close_lists()
            output.append(line.lstrip().removeprefix(RAW_LATEX_PREFIX))
            continue

        if re.match(r"^\s*<style\b", line, re.IGNORECASE):
            in_style = True
            continue
        if in_style:
            if re.search(r"</style>", line, re.IGNORECASE):
                in_style = False
            continue

        fence = re.match(r"^```\s*([A-Za-z0-9_+-]*)\s*$", line)
        if fence:
            close_lists()
            if not in_code:
                code_language = fence.group(1)
                code_lines = []
                in_code = True
            else:
                output.append(convert_code_block(code_language, code_lines))
                in_code = False
                code_language = ""
                code_lines = []
            continue

        if in_code:
            code_lines.append(line)
            continue

        html_img = re.search(r"<img\b[^>]*\bsrc=[\"']([^\"']+)[\"'][^>]*>", line, re.IGNORECASE)
        if html_img:
            close_lists()
            output.append(convert_image("", html_img.group(1)))
            continue

        if re.search(r"<score\b[^>]*></score>", line, re.IGNORECASE):
            score_image = score_event_image(current_heading)
            if score_image is not None:
                close_lists()
                output.append(score_image)
                continue

        line = clean_html_line(line)
        if not line.strip():
            close_lists()
            output.append("")
            continue

        if re.match(r"^\s*(-{3,}|\*{3,}|_{3,})\s*$", line):
            close_lists()
            continue

        heading = re.match(r"^(#{1,6})\s+(.+?)\s*$", line)
        if heading:
            close_lists()
            level = len(heading.group(1))
            command = HEADING_COMMANDS[level]
            current_heading = strip_markdown_formatting(heading.group(2))
            title = convert_inline(current_heading)
            output.append(rf"\{command}{{{title}}}")
            continue

        unordered = re.match(r"^\s*[-*]\s+(.+)$", line)
        ordered = re.match(r"^\s*\d+\.\s+(.+)$", line)
        if unordered or ordered:
            env = "itemize" if unordered else "enumerate"
            if not list_stack or list_stack[-1] != env:
                close_lists()
                output.append(rf"\begin{{{env}}}")
                list_stack.append(env)
            output.append(rf"\item {convert_inline((unordered or ordered).group(1))}")
            continue

        quote = re.match(r"^>\s?(.*)$", line)
        if quote:
            close_lists()
            output.append(r"\begin{quote}")
            output.append(convert_inline(quote.group(1)))
            output.append(r"\end{quote}")
            continue

        close_lists()
        output.append(convert_inline(line))

    close_lists()
    if in_code:
        output.append(convert_code_block(code_language, code_lines))
    return "\n".join(output)


def convert_markdown(markdown_path: Path, label: str, fallback_title: str) -> str:
    global CURRENT_DOCS_DIR, CURRENT_FOOTNOTES, CURRENT_MARKDOWN_DIR, CURRENT_MARKDOWN_FILE
    CURRENT_MARKDOWN_DIR = markdown_path.parent
    CURRENT_MARKDOWN_FILE = markdown_path
    is_index = False
    try:
        is_index = markdown_path.resolve().relative_to((CURRENT_DOCS_DIR or markdown_path.parent).resolve()) == Path("index.md")
    except ValueError:
        is_index = False
    markdown = remove_latex_false_html(strip_front_matter(markdown_path.read_text(encoding="utf-8")))
    if is_index:
        markdown = re.sub(r"^\s*!\[[^\]]*\]\([^)]*#only-(?:light|dark)\)\{[^}]*\}\s*$", "", markdown, flags=re.MULTILINE)

        def composition_card_repl(match: re.Match[str]) -> str:
            link = match.group(1)
            body = match.group(2)
            image = re.search(r"<img\b[^>]*\bsrc=[\"']([^\"']+)[\"'][^>]*>", body, re.IGNORECASE | re.DOTALL)
            if image is None:
                return ""
            figure = convert_composition_card(image.group(1), link)
            return "\n".join(RAW_LATEX_PREFIX + line for line in figure.splitlines())

        markdown = re.sub(
            r"<a\b[^>]*\bhref=[\"']([^\"']+)[\"'][^>]*>(.*?)</a>",
            composition_card_repl,
            markdown,
            flags=re.IGNORECASE | re.DOTALL,
        )
    markdown = re.sub(r"<hr\s*/?>", "---", markdown, flags=re.IGNORECASE)
    markdown = re.sub(
        r"<img\b[^>]*\bsrc=[\"']([^\"']+)[\"'][^>]*>",
        lambda match: f"![]({match.group(1)})",
        markdown,
        flags=re.IGNORECASE | re.DOTALL,
    )
    lines = markdown.splitlines()
    try:
        current_doc = markdown_path.resolve().relative_to((CURRENT_DOCS_DIR or markdown_path.parent).resolve())
    except ValueError:
        current_doc = markdown_path.name
    if current_doc == Path("integrations/index.md"):
        filtered_lines: list[str] = []
        skip_download = False
        for line in lines:
            heading = re.match(r"^(#{1,6})\s+(.+?)\s*$", line)
            if heading and strip_markdown_formatting(heading.group(2)).lower() == "download":
                skip_download = True
                continue
            if skip_download and heading and len(heading.group(1)) <= 2:
                skip_download = False
            if not skip_download:
                filtered_lines.append(line)
        lines = filtered_lines
    footnotes: dict[str, str] = {}
    content_lines: list[str] = []
    index = 0
    while index < len(lines):
        match = re.match(r"^\[\^([^\]]+)\]:\s*(.*)$", lines[index])
        if match:
            key = match.group(1)
            body = [match.group(2)]
            index += 1
            while index < len(lines) and (lines[index].startswith("    ") or lines[index].startswith("\t")):
                body.append(lines[index].strip())
                index += 1
            footnotes[key] = " ".join(part for part in body if part).strip()
            continue
        content_lines.append(lines[index])
        index += 1

    previous_footnotes = CURRENT_FOOTNOTES
    CURRENT_FOOTNOTES = footnotes
    latex = convert_markdown_lines(content_lines).strip()
    CURRENT_FOOTNOTES = previous_footnotes
    if re.search(r"^\\chapter\{", latex, flags=re.MULTILINE):
        latex = re.sub(
            r"^(\\chapter\{[^}]+\})",
            lambda match: f"{match.group(1)}\n\\label{{{label}}}",
            latex,
            count=1,
            flags=re.MULTILINE,
        )
    else:
        latex = f"\\chapter{{{latex_escape(fallback_title)}}}\n\\label{{{label}}}\n\n{latex}"
    return latex + "\n"


def write_summary(entries: list[NavEntry], output_path: Path) -> None:
    lines = [r"\chapter*{Summary}", r"\addcontentsline{toc}{chapter}{Summary}", r"\begin{itemize}"]
    current_section: tuple[str, ...] | None = None
    chapter_number = 0

    for entry in entries:
        if entry.section_path != current_section:
            if current_section:
                lines.append(r"\end{itemize}")
            current_section = entry.section_path
            if current_section:
                lines.append(rf"\item \textbf{{{latex_escape(' / '.join(current_section))}}}")
                lines.append(r"\begin{itemize}")

        if isinstance(entry, Page):
            chapter_number += 1
            label = f"chap:{slug_for_path(entry.source.relative_to(CURRENT_DOCS_DIR or entry.source.parent))}"
            lines.append(rf"\item \hyperref[{label}]{{Chapter {chapter_number} - {latex_escape(entry.title)}}}")
        else:
            lines.append(rf"\item \href{{{entry.url}}}{{{latex_escape(entry.title)}}}")

    if current_section:
            lines.append(r"\end{itemize}")
    lines.append(r"\end{itemize}")

    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def title_logo(project: dict) -> str:
    if CURRENT_DOCS_DIR is None:
        return ""
    logo = ""
    index_path = CURRENT_DOCS_DIR / "index.md"
    if index_path.exists():
        match = re.search(r"!\[[^\]]*\]\(([^)#]+)#only-light\)", index_path.read_text(encoding="utf-8"))
        if match:
            logo = match.group(1)
    if not logo:
        logo = project.get("extra", {}).get("logo_light_mode", "assets/logo.svg")
    logo_path = (CURRENT_DOCS_DIR / logo).resolve()
    if not logo_path.exists():
        return ""
    if logo_path.suffix.lower() == ".svg":
        rendered = render_svg(logo_path)
        if rendered is None:
            return ""
        logo_target = rendered
    else:
        logo_target = path_for_latex(logo_path)
    return rf"\includegraphics[width=0.22\linewidth]{{{latex_escape(logo_target)}}}\par"


def write_main(project: dict, entries: list[NavEntry], output_dir: Path) -> None:
    site_name = project.get("site_name", "Documentation")
    subtitle = "A Machine Listening System for Contemporary Music"
    copyright_text = re.sub(r"&copy;?", r"\\copyright{}", project.get("copyright", ""))
    logo = title_logo(project)
    lines = [
        r"\documentclass[11pt,oneside]{book}",
        r"\usepackage[utf8]{inputenc}",
        r"\usepackage[T1]{fontenc}",
        r"\usepackage{libertinus}",
        r"\usepackage[scaled=0.92]{inconsolata}",
        r"\usepackage{microtype}",
        r"\usepackage{geometry}",
        r"\usepackage{graphicx}",
        r"\usepackage{float}",
        r"\usepackage{xcolor}",
        r"\usepackage{amsmath}",
        r"\usepackage{amssymb}",
        r"\usepackage{hyperref}",
        r"\hypersetup{colorlinks=true, linkcolor=blue!55!black, urlcolor=blue!55!black, citecolor=blue!55!black}",
        r"\usepackage{tabularx}",
        r"\usepackage{listings}",
        r"\usepackage[most]{tcolorbox}",
        r"\tcbuselibrary{listings,breakable}",
        r"\geometry{margin=1in}",
        r"\lstset{basicstyle=\ttfamily\small,breaklines=true,columns=fullflexible,keepspaces=true}",
        r"\begin{document}",
        r"\begin{titlepage}",
        r"\centering",
        r"\vspace*{0.15\textheight}",
    ]
    if logo:
        lines.extend([logo, r"\vspace{2em}"])
    lines.extend(
        [
            rf"{{\Huge\bfseries {latex_escape(site_name)}\par}}",
            r"\vspace{1em}",
            rf"{{{latex_escape(subtitle)}\par}}",
            r"\vfill",
        ]
    )
    if copyright_text:
        lines.append(rf"{{{copyright_text}\par}}")
    lines.extend([r"\end{titlepage}", ""])
    lines.append(r"\tableofcontents")

    current_section: tuple[str, ...] = ()
    for entry in entries:
        if isinstance(entry, ExternalLink):
            continue
        if entry.section_path != current_section:
            current_section = entry.section_path
            if current_section:
                lines.append(rf"\part{{{latex_escape(' / '.join(current_section))}}}")
        lines.append(rf"\input{{chapters/{entry.tex_name}}}")

    lines.append(r"\end{document}")
    (output_dir / "main.tex").write_text("\n".join(lines) + "\n", encoding="utf-8")


def copy_assets(docs_dir: Path, output_dir: Path) -> None:
    assets = docs_dir / "assets"
    if assets.exists():
        target = output_dir / "assets"
        if target.exists():
            shutil.rmtree(target)
        shutil.copytree(assets, target)


def build_latex_project(config: Path, output_dir: Path) -> None:
    global CURRENT_CHAPTER_LINKS, CURRENT_DOCS_DIR, CURRENT_OUTPUT_DIR
    project = load_project(config)
    docs_dir = (config.parent / project.get("docs_dir", "docs")).resolve()
    CURRENT_DOCS_DIR = docs_dir
    CURRENT_OUTPUT_DIR = output_dir
    entries = flatten_nav(project.get("nav", []), docs_dir)
    pages = [entry for entry in entries if isinstance(entry, Page)]
    CURRENT_CHAPTER_LINKS = {
        page.source.resolve(): (f"chap:{slug_for_path(page.source.relative_to(docs_dir))}", page.title)
        for page in pages
    }

    chapters_dir = output_dir / "chapters"
    chapters_dir.mkdir(parents=True, exist_ok=True)
    for old_chapter in chapters_dir.glob("*.tex"):
        old_chapter.unlink()

    for page in pages:
        chapter_path = chapters_dir / page.tex_name
        label = f"chap:{slug_for_path(page.source.relative_to(docs_dir))}"
        chapter_path.write_text(convert_markdown(page.source, label, page.title), encoding="utf-8")

    summary_path = output_dir / "summary.tex"
    if summary_path.exists():
        summary_path.unlink()
    write_main(project, entries, output_dir)
    copy_assets(docs_dir, output_dir)


def compile_latex_project(output_dir: Path, runs: int = 2) -> bool:
    if shutil.which("pdflatex") is None:
        print("Skipping PDF compile: pdflatex was not found.")
        return False

    main_tex = output_dir / "main.tex"
    if not main_tex.exists():
        print(f"Skipping PDF compile: {main_tex} does not exist.")
        return False

    for pattern in ("*.aux", "*.toc", "*.out", "*.log"):
        for aux_file in output_dir.glob(pattern):
            aux_file.unlink()

    for run in range(1, runs + 1):
        result = subprocess.run(
            ["pdflatex", "-interaction=nonstopmode", "-halt-on-error", "main.tex"],
            cwd=output_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            log_path = output_dir / "main.log"
            print(f"PDF compile failed on pass {run}. See {log_path}.")
            return False

    print(f"PDF file: {output_dir / 'main.pdf'}")
    return True


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", nargs="?", default="zensical.toml", type=Path)
    parser.add_argument("-o", "--output", default=Path("build/latex"), type=Path)
    parser.add_argument("--compile", action="store_true", help="Run pdflatex after writing the LaTeX project.")
    parser.add_argument("--latex-runs", default=3, type=int, help="Number of pdflatex passes to run when compiling.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    build_latex_project(args.config.resolve(), args.output.resolve())
    print(f"Wrote LaTeX project to {args.output}")
    print(f"Main file: {args.output / 'main.tex'}")
    if args.compile:
        compile_latex_project(args.output.resolve(), max(1, args.latex_runs))


if __name__ == "__main__":
    main()
