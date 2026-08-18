#!/usr/bin/env python3
"""Convert a Zensical/MkDocs-style project into a self-contained Typst project."""

from __future__ import annotations

import argparse
import hashlib
import html
import os
import re
import shutil
import subprocess
import tomllib
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from urllib.parse import urlparse


RAW_TYPST_PREFIX = "@@RAW_TYPST_BLOCK@@"
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
CURRENT_TYPST_DIR: Path | None = None
CURRENT_CHAPTER_LINKS: dict[Path, tuple[str, str]] = {}
CURRENT_FOOTNOTES: dict[str, str] = {}
CURRENT_HEADING_BASE = 1
GRID_CARD_COUNTER = 0


@dataclass(frozen=True)
class Page:
    source: Path
    title: str
    typ_name: str
    section_path: tuple[str, ...]


@dataclass(frozen=True)
class ExternalLink:
    title: str
    url: str
    section_path: tuple[str, ...]


NavEntry = Page | ExternalLink


def typst_string(text: str) -> str:
    """Escape text for a Typst quoted string."""
    return (
        text.replace("\\", "\\\\")
        .replace('"', '\\"')
        .replace("\r", "")
        .replace("\n", "\\n")
    )


def typst_escape(text: str) -> str:
    """Escape plain text for Typst markup mode."""
    replacements = {
        "\\": r"\\",
        "#": r"\#",
        "*": r"\*",
        "_": r"\_",
        "[": r"\[",
        "]": r"\]",
        "<": r"\<",
        ">": r"\>",
        "@": r"\@",
        "$": r"\$",
        "`": r"\`",
    }
    return "".join(replacements.get(char, char) for char in text)


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
            entries.append(Page(source, title_from_markdown(source), f"{slug_for_path(Path(item))}.typ", section_path))
        elif isinstance(item, dict):
            for title, children in item.items():
                if isinstance(children, str):
                    if is_url(children):
                        entries.append(ExternalLink(title, children, section_path))
                    else:
                        source = docs_dir / children
                        entries.append(
                            Page(
                                source,
                                title_from_markdown(source),
                                f"{slug_for_path(Path(children))}.typ",
                                section_path + (title,),
                            )
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


def remove_typst_false_html(markdown: str) -> str:
    """Remove HTML blocks marked typst=false.

    For backwards compatibility, latex=false is also treated as excluded because
    existing documentation may already use that attribute to mean "web only".
    """
    false_attr = r"\b(?:typst|latex)\s*=\s*(?:\"false\"|'false'|false)(?=\s|/|>)"
    tag_name = r"[A-Za-z][A-Za-z0-9:-]*"
    block_pattern = re.compile(
        rf"<(?P<tag>{tag_name})\b(?=[^>]*{false_attr})[^>]*>.*?</(?P=tag)\s*>",
        flags=re.IGNORECASE | re.DOTALL,
    )
    previous = None
    while previous != markdown:
        previous = markdown
        markdown = block_pattern.sub("", markdown)
    return re.sub(rf"<{tag_name}\b(?=[^>]*{false_attr})[^>]*?/?>", "", markdown, flags=re.IGNORECASE)


def output_relative(path: Path) -> str:
    """Return a path relative to the Typst file currently being generated.

    Typst resolves image/include paths relative to the file containing the
    expression. Chapter files live in output/chapters, while main.typ lives
    directly in output, so using output_dir unconditionally breaks assets
    referenced from included chapter files.
    """
    base = CURRENT_TYPST_DIR or CURRENT_OUTPUT_DIR
    if base is None:
        return path.as_posix()
    return Path(os.path.relpath(path, base)).as_posix()


def asset_destination(source: Path) -> Path:
    assert CURRENT_OUTPUT_DIR is not None
    generated = CURRENT_OUTPUT_DIR / "generated-assets"
    generated.mkdir(parents=True, exist_ok=True)

    try:
        if CURRENT_DOCS_DIR is not None:
            rel = source.resolve().relative_to(CURRENT_DOCS_DIR.resolve())
            return generated / rel
    except ValueError:
        pass

    digest = hashlib.sha256(str(source.resolve()).encode("utf-8")).hexdigest()[:16]
    return generated / f"{digest}-{source.name}"


def copy_local_asset(source: Path) -> str | None:
    if CURRENT_OUTPUT_DIR is None or not source.exists() or not source.is_file():
        return None
    target = asset_destination(source)
    target.parent.mkdir(parents=True, exist_ok=True)
    if not target.exists() or source.stat().st_mtime_ns > target.stat().st_mtime_ns:
        shutil.copy2(source, target)
    return output_relative(target)


def download_remote_image(url: str) -> str | None:
    if CURRENT_OUTPUT_DIR is None:
        return None
    suffix = Path(urlparse(url).path).suffix.lower()
    if suffix not in {".png", ".jpg", ".jpeg", ".svg", ".gif", ".webp", ".pdf"}:
        return None
    generated_dir = CURRENT_OUTPUT_DIR / "generated-images" / "remote"
    generated_dir.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256(url.encode("utf-8")).hexdigest()[:16]
    image_path = generated_dir / f"{digest}{suffix}"
    if not image_path.exists():
        try:
            with urllib.request.urlopen(url, timeout=30) as response:
                image_path.write_bytes(response.read())
        except Exception as exc:
            print(f"Warning: could not download {url}: {exc}")
            return None
    return output_relative(image_path)


def prepare_mermaid_source(source: str) -> str:
    lines = [line for line in source.splitlines() if not re.match(r"^\s*%%\{init:.*\}%%\s*$", line)]
    return "\n".join([MERMAID_INIT, *lines]).strip()


def render_mermaid(source: str) -> str | None:
    if CURRENT_OUTPUT_DIR is None or shutil.which("mmdc") is None:
        return None
    generated_dir = CURRENT_OUTPUT_DIR / "generated-images" / "mermaid"
    generated_dir.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256(source.encode("utf-8")).hexdigest()[:16]
    mermaid_path = generated_dir / f"{digest}.mmd"
    png_path = generated_dir / f"{digest}.png"
    if png_path.exists():
        return output_relative(png_path)

    mermaid_path.write_text(source, encoding="utf-8")
    command = ["mmdc", "-i", str(mermaid_path), "-o", str(png_path), "-b", "white", "-s", "2"]

    chromium = shutil.which("chromium") or shutil.which("chromium-browser") or shutil.which("google-chrome")
    if chromium:
        puppeteer_config = generated_dir / "puppeteer-config.json"
        puppeteer_config.write_text(
            '{"executablePath": "' + chromium.replace("\\", "\\\\") + '", "args": ["--no-sandbox", "--disable-setuid-sandbox"]}\n',
            encoding="utf-8",
        )
        command.extend(["-p", str(puppeteer_config)])

    result = subprocess.run(
        command,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    if result.returncode != 0 or not png_path.exists():
        return None
    return output_relative(png_path)


def raw_block(source: str, language: str = "") -> str:
    lang = f', lang: "{typst_string(language)}"' if language else ""
    return f'#raw("{typst_string(source)}", block: true{lang})'


def convert_mermaid(source_lines: list[str]) -> str:
    source = prepare_mermaid_source("\n".join(source_lines).strip())
    rendered = render_mermaid(source)
    if rendered is not None:
        return f'#align(center)[#image("{typst_string(rendered)}", width: 90%)]'
    return raw_block("\n".join(source_lines), "mermaid")


def convert_code_block(language: str, source_lines: list[str]) -> str:
    if language.lower() == "mermaid":
        return convert_mermaid(source_lines)
    return raw_block("\n".join(source_lines), language)


def resolve_local_markdown_link(target: str) -> tuple[str, str] | None:
    if CURRENT_MARKDOWN_DIR is None or CURRENT_MARKDOWN_FILE is None:
        return None
    target = target.split("{", 1)[0].strip()
    target_path = target.split("#", 1)[0]
    if not target_path:
        return CURRENT_CHAPTER_LINKS.get(CURRENT_MARKDOWN_FILE.resolve())

    candidate = (CURRENT_MARKDOWN_DIR / target_path).resolve()
    candidates = [candidate]

    suffix = Path(target_path).suffix.lower()
    if suffix == "":
        candidates.extend([candidate.with_suffix(".md"), candidate / "index.md"])
    elif suffix != ".md":
        return None

    for resolved in candidates:
        chapter = CURRENT_CHAPTER_LINKS.get(resolved)
        if chapter is not None:
            return chapter
    return None


def resolve_image_target(target: str) -> str | None:
    target = target.split("#", 1)[0].strip()
    if not target:
        return None
    if is_url(target):
        return download_remote_image(target)
    if CURRENT_MARKDOWN_DIR is None:
        return target
    resolved = (CURRENT_MARKDOWN_DIR / target).resolve()
    copied = copy_local_asset(resolved)
    return copied if copied is not None else target.lstrip("./")


def convert_image(alt: str, target: str) -> str:
    resolved = resolve_image_target(target)
    if not resolved:
        return ""
    caption = strip_markdown_formatting(alt)
    image_expr = f'image("{typst_string(resolved)}", width: 80%)'
    if caption:
        return f'#figure({image_expr}, caption: [{convert_inline(caption)}])'
    return f'#align(center)[#{image_expr}]'


def convert_card_image(target: str) -> str:
    resolved = resolve_image_target(target)
    if not resolved:
        return ""
    return f'#align(center)[#image("{typst_string(resolved)}", width: 100%, height: 72mm, fit: "contain")]'


def convert_composition_card(image_target: str, link_target: str) -> str:
    resolved = resolve_image_target(image_target)
    if not resolved:
        return ""
    return (
        '#figure(\n'
        f'  image("{typst_string(resolved)}", width: 80%),\n'
        f'  caption: [#link("{typst_string(link_target)}")[{typst_escape(link_target)}]],\n'
        ')'
    )


def latex_math_to_typst(source: str) -> str:
    """Best-effort conversion from common LaTeX math to Typst math."""
    text = source.strip()

    # Convert brace-delimited commands before generic substitutions. A regex like
    # ``\\frac\{([^{}]+)\}\{([^{}]+)\}`` fails as soon as either argument contains
    # nested braces (for example ``S_{n,k}``).
    def read_group(value: str, start: int) -> tuple[str, int] | None:
        if start >= len(value) or value[start] != "{":
            return None
        depth = 0
        for i in range(start, len(value)):
            if value[i] == "{":
                depth += 1
            elif value[i] == "}":
                depth -= 1
                if depth == 0:
                    return value[start + 1:i], i + 1
        return None

    def replace_two_group_command(value: str, command: str, build) -> str:
        needle = "\\" + command
        out: list[str] = []
        i = 0
        while i < len(value):
            if value.startswith(needle, i):
                j = i + len(needle)
                while j < len(value) and value[j].isspace():
                    j += 1
                first = read_group(value, j)
                if first is not None:
                    a, j2 = first
                    while j2 < len(value) and value[j2].isspace():
                        j2 += 1
                    second = read_group(value, j2)
                    if second is not None:
                        b, j3 = second
                        out.append(build(latex_math_to_typst(a), latex_math_to_typst(b)))
                        i = j3
                        continue
            out.append(value[i])
            i += 1
        return "".join(out)

    def replace_one_group_command(value: str, command: str, build) -> str:
        needle = "\\" + command
        out: list[str] = []
        i = 0
        while i < len(value):
            if value.startswith(needle, i):
                j = i + len(needle)
                while j < len(value) and value[j].isspace():
                    j += 1
                group = read_group(value, j)
                if group is not None:
                    body, j2 = group
                    out.append(build(latex_math_to_typst(body)))
                    i = j2
                    continue
            out.append(value[i])
            i += 1
        return "".join(out)

    # Repeat so nested fractions are handled from the outside through recursion.
    text = replace_two_group_command(text, "frac", lambda a, b: f"({a})/({b})")
    text = replace_one_group_command(text, "text", lambda a: f'"{a}"')
    text = replace_one_group_command(text, "mathrm", lambda a: f"upright({a})")
    text = replace_one_group_command(text, "mathbf", lambda a: f"bold({a})")

    replacements = {
        r"\alpha": "alpha", r"\beta": "beta", r"\gamma": "gamma", r"\delta": "delta",
        r"\epsilon": "epsilon", r"\varepsilon": "epsilon.alt", r"\zeta": "zeta",
        r"\eta": "eta", r"\theta": "theta", r"\vartheta": "theta.alt", r"\iota": "iota",
        r"\kappa": "kappa", r"\lambda": "lambda", r"\mu": "mu", r"\nu": "nu",
        r"\xi": "xi", r"\pi": "pi", r"\rho": "rho", r"\sigma": "sigma",
        r"\tau": "tau", r"\upsilon": "upsilon", r"\phi": "phi", r"\varphi": "phi.alt",
        r"\chi": "chi", r"\psi": "psi", r"\omega": "omega",
        r"\Gamma": "Gamma", r"\Delta": "Delta", r"\Theta": "Theta", r"\Lambda": "Lambda",
        r"\Xi": "Xi", r"\Pi": "Pi", r"\Sigma": "Sigma", r"\Phi": "Phi",
        r"\Psi": "Psi", r"\Omega": "Omega",
        r"\times": "times", r"\cdot": "dot", r"\pm": "+-", r"\mp": "-+",
        r"\leq": "<=", r"\geq": ">=", r"\neq": "!=", r"\approx": "approx",
        r"\infty": "infinity", r"\sum": "sum", r"\prod": "product", r"\int": "integral",
        r"\partial": "diff", r"\rightarrow": "->", r"\leftarrow": "<-",
        r"\Rightarrow": "=>", r"\Leftarrow": "<=", r"\in": "in", r"\notin": "in.not",
        r"\subset": "subset", r"\subseteq": "subset.eq", r"\cup": "union", r"\cap": "sect",
        r"\ldots": "dots.h", r"\cdots": "dots.h.c", r"\sqrt": "sqrt",
        # TeX function commands must lose their backslash. Leaving ``\\log`` in
        # Typst produces a stray backslash and can lead to errors such as
        # ``unknown variable: og``.
        r"\log": "log", r"\ln": "ln", r"\exp": "exp", r"\lim": "lim",
        r"\sin": "sin", r"\cos": "cos", r"\tan": "tan", r"\cot": "cot",
        r"\sec": "sec", r"\csc": "csc", r"\sinh": "sinh", r"\cosh": "cosh",
        r"\tanh": "tanh", r"\min": "min", r"\max": "max", r"\det": "det",
    }
    # Longest commands first so e.g. \subseteq is handled before \subset.
    for old in sorted(replacements, key=len, reverse=True):
        text = text.replace(old, replacements[old])

    text = text.replace(r"\left", "").replace(r"\right", "")
    text = text.replace("{", "(").replace("}", ")")

    # Preserve quoted text while splitting TeX-style juxtaposed letters. In TeX,
    # RMS means R M S; Typst otherwise interprets RMS as one identifier.
    quoted: list[str] = []

    def stash_quoted(match: re.Match[str]) -> str:
        quoted.append(match.group(0))
        return f"@@MATHQUOTED{len(quoted) - 1}@@"

    text = re.sub(r'"(?:\\.|[^"\\])*"', stash_quoted, text)

    known_math_names = {
        "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "eta",
        "theta", "iota", "kappa", "lambda", "mu", "nu", "xi", "pi",
        "rho", "sigma", "tau", "upsilon", "phi", "chi", "psi", "omega",
        "Gamma", "Delta", "Theta", "Lambda", "Xi", "Pi", "Sigma", "Phi",
        "Psi", "Omega", "times", "dot", "approx", "infinity", "sum",
        "product", "integral", "diff", "in", "subset", "union", "sect",
        "sqrt", "upright", "italic", "bold", "sin", "cos", "tan", "cot",
        "sec", "csc", "sinh", "cosh", "tanh", "log", "ln", "exp", "lim",
        "min", "max", "mod", "gcd", "lcm", "det", "Pr",
    }

    def split_identifier(match: re.Match[str]) -> str:
        word = match.group(0)
        if word in known_math_names or word.startswith("MATHQUOTED"):
            return word
        if len(word) <= 1:
            return word
        return " ".join(word)

    text = re.sub(r"(?<![@.A-Za-z0-9])[A-Za-z]+(?![A-Za-z0-9])", split_identifier, text)

    for index, value in enumerate(quoted):
        text = text.replace(f"@@MATHQUOTED{index}@@", value)

    return text

def convert_math_token(token: str) -> str:
    if token.startswith("$$") and token.endswith("$$"):
        body = token[2:-2]
        return f'$ {latex_math_to_typst(body)} $'
    if token.startswith(r"\[") and token.endswith(r"\]"):
        body = token[2:-2]
        return f'$ {latex_math_to_typst(body)} $'
    if token.startswith(r"\(") and token.endswith(r"\)"):
        body = token[2:-2]
        return f'${latex_math_to_typst(body)}$'
    if token.startswith("$") and token.endswith("$"):
        body = token[1:-1]
        return f'${latex_math_to_typst(body)}$'
    return typst_escape(token)


def convert_inline(text: str) -> str:
    placeholders: list[str] = []

    def stash(value: str) -> str:
        placeholders.append(value)
        return f"@@PLACEHOLDER{len(placeholders) - 1}@@"

    def stash_math(match: re.Match[str]) -> str:
        return stash(convert_math_token(match.group(0)))

    text = html.unescape(text)
    text = re.sub(r"\$\$.*?\$\$", stash_math, text, flags=re.DOTALL)
    text = re.sub(r"\\\[.*?\\\]", stash_math, text, flags=re.DOTALL)
    text = re.sub(r"\\\(.*?\\\)", stash_math, text, flags=re.DOTALL)
    text = re.sub(r"(?<!\\)\$(?!\s)(.+?)(?<!\s)(?<!\\)\$", stash_math, text)

    def code_repl(match: re.Match[str]) -> str:
        return stash(f'#raw("{typst_string(match.group(1))}")')

    def bold_repl(match: re.Match[str]) -> str:
        return stash(f'*{convert_inline(match.group(1))}*')

    def italic_repl(match: re.Match[str]) -> str:
        return stash(f'_{convert_inline(match.group(1))}_')

    def image_repl(match: re.Match[str]) -> str:
        return stash(convert_image(match.group(1), match.group(2)))

    def footnote_repl(match: re.Match[str]) -> str:
        footnote = CURRENT_FOOTNOTES.get(match.group(1))
        if footnote is None:
            return typst_escape(match.group(0))
        return stash(f'#footnote[{convert_inline(footnote)}]')

    def link_repl(match: re.Match[str]) -> str:
        label = convert_inline(match.group(1))
        target = match.group(2).split("{", 1)[0].strip()
        if is_url(target):
            return stash(f'#link("{typst_string(target)}")[{label}]')
        local_chapter = resolve_local_markdown_link(target)
        if local_chapter is not None:
            chapter_label, chapter_title = local_chapter
            if not label:
                label = typst_escape(chapter_title)
            return stash(f'#link(<{chapter_label}>)[{label}]')
        suffix = f' (#raw("{typst_string(target)}"))' if target else ""
        return stash(label + suffix)

    text = re.sub(r"!\[([^\]]*)\]\(([^)]+)\)(?:\{[^}]*\})?", image_repl, text)
    text = re.sub(r"\[([^\]]*)\]\(([^)]+)\)(?:\{[^}]*\})?", link_repl, text)
    text = re.sub(r"\[\^([^\]]+)\]", footnote_repl, text)
    text = re.sub(r"`([^`]+)`", code_repl, text)
    text = re.sub(r"\*\*(.+?)\*\*", bold_repl, text)
    text = re.sub(r"__(.+?)__", bold_repl, text)
    text = re.sub(r"(?<!\*)\*([^*\n]+)\*(?!\*)", italic_repl, text)
    text = re.sub(r"(?<!_)_([^_\n]+)_(?!_)", italic_repl, text)
    text = re.sub(r":[-+a-zA-Z0-9_]+:", "", text)

    escaped = typst_escape(text)
    changed = True
    while changed:
        changed = False
        for index, value in enumerate(placeholders):
            raw_placeholder = f"@@PLACEHOLDER{index}@@"
            escaped_placeholder = typst_escape(raw_placeholder)
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
    line = re.sub(r"<br\s*/?>", "  ", line, flags=re.IGNORECASE)
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
    copied = copy_local_asset(image_path)
    if copied is None:
        return None
    return f'#align(center)[#image("{typst_string(copied)}", width: 80%)]'


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


def match_code_fence(line: str) -> re.Match[str] | None:
    return re.match(r"^```\s*([A-Za-z0-9_+-]*)[^\n`]*$", line)


def convert_admonition(kind: str, title: str, body: list[str]) -> str:
    heading = title or kind.title()
    converted_body = convert_markdown_lines(body).strip()
    if not converted_body:
        converted_body = "#v(1mm)"
    return RAW_TYPST_PREFIX + "\n".join(
        [
            "#block(",
            "  width: 100%,",
            "  inset: 8pt,",
            "  radius: 4pt,",
            "  stroke: 0.5pt + luma(70%),",
            "  fill: luma(97%),",
            ")[",
            f"  *{convert_inline(heading)}*",
            "",
            converted_body,
            "]",
            "",
        ]
    )


def strip_card_body_indent(lines: list[str]) -> list[str]:
    stripped = [line[4:] if line.startswith("    ") else line[1:] if line.startswith("\t") else line for line in lines]
    while stripped and not stripped[0].strip():
        stripped.pop(0)
    while stripped and not stripped[-1].strip():
        stripped.pop()
    return stripped


def parse_grid_card_items(lines: list[str]) -> list[tuple[str, list[str]]]:
    cards: list[tuple[str, list[str]]] = []
    title: str | None = None
    body: list[str] = []
    in_code = False

    def flush() -> None:
        nonlocal title, body
        if title is not None:
            cards.append((title, strip_card_body_indent(body)))
        title = None
        body = []

    for line in lines:
        if re.match(r"^\s*```\s*", line):
            in_code = not in_code
        item = re.match(r"^-\s+(.+?)\s*$", line)
        if item and not in_code:
            flush()
            title = strip_markdown_formatting(item.group(1))
            continue
        if title is not None:
            body.append(line)
    flush()
    return cards


def convert_card_body(lines: list[str]) -> str:
    output: list[str] = []
    index = 0
    while index < len(lines):
        line = lines[index].rstrip()
        if not line.strip():
            index += 1
            continue

        fence = match_code_fence(line)
        if fence:
            language = fence.group(1)
            code_lines: list[str] = []
            index += 1
            while index < len(lines) and not match_code_fence(lines[index]):
                code_lines.append(lines[index].rstrip())
                index += 1
            if index < len(lines):
                index += 1
            output.append(convert_code_block(language, [sanitize_code_line(x) for x in code_lines]))
            continue

        image = re.match(r"^!\[[^\]]*\]\(([^)]+)\)(?:\{[^}]*\})?\s*$", line)
        if image:
            output.append(convert_card_image(image.group(1)))
            index += 1
            continue

        paragraph = [line.strip()]
        index += 1
        while (
            index < len(lines)
            and lines[index].strip()
            and not match_code_fence(lines[index])
            and not re.match(r"^!\[[^\]]*\]\(([^)]+)\)", lines[index].strip())
        ):
            paragraph.append(lines[index].strip())
            index += 1
        output.append(convert_inline(" ".join(paragraph)))

    return "\n\n".join(output) if output else "#v(1mm)"


def convert_grid_cards(block_lines: list[str]) -> str:
    global GRID_CARD_COUNTER
    cards = parse_grid_card_items(block_lines)
    if not cards:
        return ""
    GRID_CARD_COUNTER += 1

    columns = 1 if len(cards) <= 1 else 2
    column_spec = "(1fr,)" if columns == 1 else "(1fr, 1fr)"
    items: list[str] = []
    for title, body in cards:
        items.append(
            "\n".join(
                [
                    "block(",
                    "  width: 100%,",
                    "  inset: 8pt,",
                    "  radius: 4pt,",
                    "  stroke: 0.35pt + luma(82%),",
                    ")[",
                    f"*{convert_inline(title)}*",
                    "",
                    convert_card_body(body),
                    "]",
                ]
            )
        )

    return RAW_TYPST_PREFIX + "#grid(\n  columns: " + column_spec + ",\n  gutter: 8pt,\n  " + ",\n  ".join(items) + ",\n)"


def normalize_grid_cards(lines: list[str]) -> list[str]:
    output: list[str] = []
    index = 0
    in_code = False
    start_pattern = re.compile(r"<div\b(?=[^>]*\bgrid\b)(?=[^>]*\bcards\b)(?=[^>]*\bmarkdown\b)[^>]*>", re.IGNORECASE)

    while index < len(lines):
        line = lines[index]
        if match_code_fence(line):
            in_code = not in_code
            output.append(line)
            index += 1
            continue

        if not in_code and start_pattern.search(line):
            index += 1
            block: list[str] = []
            while index < len(lines) and not re.search(r"</div\s*>", lines[index], re.IGNORECASE):
                block.append(lines[index])
                index += 1
            if index < len(lines):
                index += 1
            output.append(convert_grid_cards(block))
            continue

        output.append(line)
        index += 1

    return output


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

    def normalize_cells(cells: list[str]) -> list[str]:
        padded = cells[:column_count] + [""] * max(0, column_count - len(cells))
        return padded[:column_count]

    cells: list[str] = []
    for cell in normalize_cells(header):
        cells.append(f"[*{convert_inline(cell)}*]")
    for row in body_rows:
        for cell in normalize_cells(row):
            cells.append(f"[{convert_inline(cell)}]")

    return RAW_TYPST_PREFIX + "\n".join(
        [
            "#table(",
            f"  columns: {column_count},",
            "  stroke: 0.4pt + luma(75%),",
            "  inset: 5pt,",
            "  align: left + top,",
            "  " + ",\n  ".join(cells) + ",",
            ")",
            "",
        ]
    )


def normalize_tables(lines: list[str]) -> list[str]:
    output: list[str] = []
    index = 0
    in_code = False

    while index < len(lines):
        line = lines[index]
        if match_code_fence(line):
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
        match = re.match(r'^([!?]{3})\s+([A-Za-z0-9_-]+)(?:\s+"([^"]+)")?\s*$', lines[index])
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


def heading_markup(level: int, title: str) -> str:
    effective = max(1, min(6, CURRENT_HEADING_BASE + level - 1))
    return "=" * effective + " " + title


def convert_markdown_lines(lines: list[str]) -> str:
    lines = normalize_grid_cards(lines)
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
        if list_stack:
            output.append("")
            list_stack.clear()

    for raw_line in lines:
        line = raw_line.rstrip()

        if line.lstrip().startswith(RAW_TYPST_PREFIX):
            close_lists()
            output.append(line.lstrip().removeprefix(RAW_TYPST_PREFIX))
            continue

        if re.match(r"^\s*<style\b", line, re.IGNORECASE):
            in_style = True
            continue
        if in_style:
            if re.search(r"</style>", line, re.IGNORECASE):
                in_style = False
            continue

        fence = match_code_fence(line)
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
            code_lines.append(sanitize_code_line(line))
            continue

        html_img = re.search(r'<img\b[^>]*\bsrc=["\']([^"\']+)["\'][^>]*>', line, re.IGNORECASE)
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
            output.append("#line(length: 100%, stroke: 0.5pt + luma(75%))")
            continue

        heading = re.match(r"^(#{1,6})\s+(.+?)\s*$", line)
        if heading:
            close_lists()
            level = len(heading.group(1))
            current_heading = strip_markdown_formatting(heading.group(2))
            output.append(heading_markup(level, convert_inline(current_heading)))
            continue

        unordered = re.match(r"^\s*[-*]\s+(.+)$", line)
        ordered = re.match(r"^\s*\d+\.\s+(.+)$", line)
        if unordered or ordered:
            env = "unordered" if unordered else "ordered"
            if list_stack and list_stack[-1] != env:
                close_lists()
            if not list_stack:
                list_stack.append(env)
            marker = "-" if env == "unordered" else "+"
            output.append(f"{marker} {convert_inline((unordered or ordered).group(1))}")
            continue

        quote = re.match(r"^>\s?(.*)$", line)
        if quote:
            close_lists()
            output.append(
                "#block(inset: (left: 10pt), stroke: (left: 2pt + luma(65%)))["
                + convert_inline(quote.group(1))
                + "]"
            )
            continue

        close_lists()
        output.append(convert_inline(line))

    close_lists()
    if in_code:
        output.append(convert_code_block(code_language, code_lines))
    return "\n".join(output)


def convert_markdown(markdown_path: Path, label: str, fallback_title: str, heading_base: int = 1) -> str:
    global CURRENT_DOCS_DIR, CURRENT_FOOTNOTES, CURRENT_HEADING_BASE, CURRENT_MARKDOWN_DIR, CURRENT_MARKDOWN_FILE
    CURRENT_MARKDOWN_DIR = markdown_path.parent
    CURRENT_MARKDOWN_FILE = markdown_path
    CURRENT_HEADING_BASE = heading_base

    is_index = False
    try:
        is_index = markdown_path.resolve().relative_to((CURRENT_DOCS_DIR or markdown_path.parent).resolve()) == Path("index.md")
    except ValueError:
        is_index = False

    markdown = remove_typst_false_html(strip_front_matter(markdown_path.read_text(encoding="utf-8")))
    if is_index:
        markdown = re.sub(
            r"^\s*!\[[^\]]*\]\([^)]*#only-(?:light|dark)\)\{[^}]*\}\s*$",
            "",
            markdown,
            flags=re.MULTILINE,
        )

        def composition_card_repl(match: re.Match[str]) -> str:
            link = match.group(1)
            body = match.group(2)
            image = re.search(r'<img\b[^>]*\bsrc=["\']([^"\']+)["\'][^>]*>', body, re.IGNORECASE | re.DOTALL)
            if image is None:
                return ""
            figure = convert_composition_card(image.group(1), link)
            return "\n".join(RAW_TYPST_PREFIX + line for line in figure.splitlines())

        markdown = re.sub(
            r'<a\b[^>]*\bhref=["\']([^"\']+)["\'][^>]*>(.*?)</a>',
            composition_card_repl,
            markdown,
            flags=re.IGNORECASE | re.DOTALL,
        )

    markdown = re.sub(r"<hr\s*/?>", "---", markdown, flags=re.IGNORECASE)
    markdown = re.sub(
        r'<img\b[^>]*\bsrc=["\']([^"\']+)["\'][^>]*>',
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
    typst = convert_markdown_lines(content_lines).strip()
    CURRENT_FOOTNOTES = previous_footnotes

    heading_pattern = re.compile(r"^(=+)\s+(.+)$", flags=re.MULTILINE)
    first_heading = heading_pattern.search(typst)
    if first_heading:
        line = first_heading.group(0)
        replacement = line + f" <{label}>"
        typst = typst[: first_heading.start()] + replacement + typst[first_heading.end() :]
    else:
        typst = f'{"=" * heading_base} {typst_escape(fallback_title)} <{label}>\n\n{typst}'

    return typst + "\n"


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
    copied = copy_local_asset(logo_path)
    if copied is None:
        return ""
    return f'#image("{typst_string(copied)}", width: 22%)'


def section_heading_line(section_path: tuple[str, ...]) -> str:
    title = " / ".join(section_path)
    return f'= {typst_escape(title)}'


def write_main(project: dict, entries: list[NavEntry], output_dir: Path) -> None:
    global CURRENT_TYPST_DIR
    CURRENT_TYPST_DIR = output_dir
    site_name = project.get("site_name", "Documentation")
    subtitle = project.get("extra", {}).get("subtitle", "A Machine Listening System for Contemporary Music")
    copyright_text = html.unescape(re.sub(r"&copy;?", "©", project.get("copyright", "")))
    logo = title_logo(project)

    lines = [
        '#set page(paper: "a4", margin: 25mm)',
        '#set text(size: 11pt)',
        '#set par(justify: true, leading: 0.7em)',
        '#set heading(numbering: "1.1")',
        '#show heading.where(level: 1): it => block(above: 1.4em, below: 0.8em, text(size: 20pt, weight: "bold", it))',
        '#show heading.where(level: 2): it => block(above: 1.2em, below: 0.6em, text(size: 16pt, weight: "bold", it))',
        '#show raw.where(block: true): it => block(fill: luma(96%), inset: 7pt, radius: 3pt, width: 100%, it)',
        '',
        '#pagebreak()',
        '#align(center)[',
        '  #v(35mm)',
    ]
    if logo:
        lines.extend([f"  {logo}", "  #v(16pt)"])
    lines.extend(
        [
            f'  #text(size: 26pt, weight: "bold")[{typst_escape(site_name)}]',
            '  #v(8pt)',
            f'  #text(size: 13pt)[{typst_escape(subtitle)}]',
            '  #v(20mm)',
        ]
    )
    if copyright_text:
        lines.append(f'  #text(size: 9pt)[{typst_escape(copyright_text)}]')
    lines.extend([
        ']',
        '#pagebreak()',
        '#outline(title: [Contents])',
        '#pagebreak()',
        '',
    ])

    current_section: tuple[str, ...] = ()
    for entry in entries:
        if isinstance(entry, ExternalLink):
            continue
        if entry.section_path != current_section:
            current_section = entry.section_path
            if current_section:
                lines.extend(["#pagebreak()", section_heading_line(current_section), ""])
        lines.append(f'#include "chapters/{typst_string(entry.typ_name)}"')
        lines.append("")

    (output_dir / "main.typ").write_text("\n".join(lines) + "\n", encoding="utf-8")


def copy_assets(docs_dir: Path, output_dir: Path) -> None:
    assets = docs_dir / "assets"
    if assets.exists():
        target = output_dir / "assets"
        if target.exists():
            shutil.rmtree(target)
        shutil.copytree(assets, target)


def build_typst_project(config: Path, output_dir: Path) -> None:
    global CURRENT_CHAPTER_LINKS, CURRENT_DOCS_DIR, CURRENT_OUTPUT_DIR, CURRENT_TYPST_DIR, GRID_CARD_COUNTER
    project = load_project(config)
    docs_dir = (config.parent / project.get("docs_dir", "docs")).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    CURRENT_DOCS_DIR = docs_dir
    CURRENT_OUTPUT_DIR = output_dir
    GRID_CARD_COUNTER = 0

    entries = flatten_nav(project.get("nav", []), docs_dir)
    pages = [entry for entry in entries if isinstance(entry, Page)]
    CURRENT_CHAPTER_LINKS = {
        page.source.resolve(): (f"chap_{slug_for_path(page.source.relative_to(docs_dir))}", page.title)
        for page in pages
    }

    chapters_dir = output_dir / "chapters"
    chapters_dir.mkdir(parents=True, exist_ok=True)
    for old_chapter in chapters_dir.glob("*.typ"):
        old_chapter.unlink()

    copy_assets(docs_dir, output_dir)

    for page in pages:
        chapter_path = chapters_dir / page.typ_name
        CURRENT_TYPST_DIR = chapter_path.parent
        label = f"chap_{slug_for_path(page.source.relative_to(docs_dir))}"
        heading_base = 2 if page.section_path else 1
        chapter_path.write_text(
            convert_markdown(page.source, label, page.title, heading_base=heading_base),
            encoding="utf-8",
        )

    write_main(project, entries, output_dir)


def compile_typst_project(output_dir: Path) -> bool:
    typst = shutil.which("typst")
    if typst is None:
        print("Skipping PDF compile: typst was not found.")
        return False

    main_typ = output_dir / "main.typ"
    if not main_typ.exists():
        print(f"Skipping PDF compile: {main_typ} does not exist.")
        return False

    result = subprocess.run(
        [typst, "compile", "main.typ", "main.pdf"],
        cwd=output_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        print("Typst compilation failed:")
        print(result.stdout)
        return False

    print(f"PDF file: {output_dir / 'main.pdf'}")
    return True


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("config", nargs="?", default="zensical.toml", type=Path)
    parser.add_argument("-o", "--output", default=Path("build/typst"), type=Path)
    parser.add_argument("--compile", action="store_true", help="Run `typst compile` after writing the Typst project.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    output_dir = args.output.resolve()
    build_typst_project(args.config.resolve(), output_dir)
    print(f"Wrote Typst project to {output_dir}")
    print(f"Main file: {output_dir / 'main.typ'}")
    if args.compile:
        compile_typst_project(output_dir)


if __name__ == "__main__":
    main()
