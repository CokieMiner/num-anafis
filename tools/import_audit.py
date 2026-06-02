#!/usr/bin/env python3
"""Audit Rust imports and module-boundary usage.

Checks:
- inline qualified paths used in code bodies outside `use` statements
- deep relative imports such as `use super::super::...`
- cross-boundary imports into `logic/`
- surface bypasses that jump past a boundary's top-level `logic` re-export
- emits a Graphviz DOT graph with staircase-boundary clusters

Usage:
    python3 tools/import_audit.py
    python3 tools/import_audit.py --json
    python3 tools/import_audit.py --dot tools/import_graph.dot
"""

from __future__ import annotations

import argparse
import json
import keyword
import re
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path

# ============================================================================
# CONFIGURATION
# ============================================================================
# Project-specific exceptions can be tuned here without digging into the logic.

CONFIG = {
    # Imports from `crate::` that are globally permitted (e.g., global constants).
    "whitelisted_crate_root_imports": {
        "crate::EPSILON",
        "crate::DEFAULT_MAX_DEPTH",
        "crate::DEFAULT_MAX_NODES",
    },
    
    # Specific type aliases that are allowed without warnings.
    "whitelisted_alias_imports": {
        ("FmtResult", "std::fmt::Result"),
        ("backend", "f32_ops"),
        ("backend", "f64_ops"),
        ("backend", "rug_ops"),
        ("backend", "i32_math"),
        ("backend", "i64_math"),
    },
    
    # Module-specific allowed alias imports (e.g., Python bindings renaming things).
    "whitelisted_alias_import_prefixes": {
        "src/bindings/python/": {
            ("RustExpr", "src::core::Expr"),
            ("RustSymbol", "src::core::Symbol"),
            ("RustContext", "src::core::Context"),
            ("RustVmEvaluator", "src::evaluator::VmEvaluator"),
            ("parse_expr", "src::parser::parse"),
            ("rust_eval_f64", "src::evaluator::eval_f64"),
            ("rust_diff", "src::diff::diff"),
            ("rust_simplify", "src::simplification::simplify"),
        },
    },

    # Primitive types that act as namespaces but are not PascalCase.
    "rust_primitive_types": {
        "f32", "f64", "f64x4", "i8", "i16", "i32", "i64", "i128", "isize",
        "u8", "u16", "u32", "u64", "u128", "usize", "str", "bool", "char"
    },

    # Specific inline qualified paths that are allowed (e.g. array::from_fn)
    "whitelisted_inline_qualified_paths": {
        "array::from_fn",
    }
}
# ============================================================================

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
RUST_KEYWORDS = {
    "Self", "async", "await", "break", "const", "continue", "crate", "dyn",
    "else", "enum", "extern", "false", "fn", "for", "if", "impl", "in", "let",
    "loop", "match", "mod", "move", "mut", "pub", "ref", "return", "self",
    "static", "struct", "super", "trait", "true", "type", "unsafe", "use",
    "where", "while",
}

USE_RE = re.compile(r"^\s*(?:pub(?:\([^)]+\))?\s+)?use\s+(.+?);$")
DEEP_RELATIVE_RE = re.compile(r"^\s*(?:pub(?:\([^)]+\))?\s+)?use\s+super::super(?:::super)*::")
INLINE_MOD_RE = re.compile(r"^\s*(?:pub(?:\([^)]+\))?\s+)?mod\s+([A-Za-z_][A-Za-z0-9_]*)\s*\{")
MACRO_RULES_RE = re.compile(r"^\s*(?:#\[.*?\]\s*)*macro_rules!\s+[A-Za-z_][A-Za-z0-9_]*\s*\{")
QUALIFIED_PATH_RE = re.compile(r"(?<!\$)\b([A-Za-z_][A-Za-z0-9_]*(?:::[A-Za-z_][A-Za-z0-9_]*)+)")
MACRO_INVOCATION_RE = re.compile(r"\b([A-Za-z_][A-Za-z0-9_]*)\s*!\s*[(\{\[]")
TEST_PATH_PARTS = {"tests", "benches"}


@dataclass
class Finding:
    kind: str
    path: str
    line: int
    text: str

@dataclass
class ImportEdge:
    source: str
    target: str
    raw: str
    path: str
    line: int
    internal: bool

@dataclass
class AliasImport:
    path: str
    line: int
    source_module: str
    target: str
    alias: str
    usage_count: int


def clean_rust_code(text: str) -> str:
    result = []
    i = 0
    n = len(text)
    state = "NORMAL"
    block_comment_depth = 0
    raw_string_hashes = 0
    
    while i < n:
        if state == "NORMAL":
            if text.startswith("//", i):
                state = "LINE_COMMENT"
                result.append("  ")
                i += 2
            elif text.startswith("/*", i):
                state = "BLOCK_COMMENT"
                block_comment_depth = 1
                result.append("  ")
                i += 2
            elif text.startswith('r#', i) or text.startswith('r"', i):
                hashes = 0
                idx = i + 1
                while idx < n and text[idx] == '#':
                    hashes += 1
                    idx += 1
                if idx < n and text[idx] == '"':
                    state = "RAW_STRING"
                    raw_string_hashes = hashes
                    result.append(" " * (idx - i + 1))
                    i = idx + 1
                else:
                    result.append(text[i])
                    i += 1
            elif text[i] == '"':
                state = "STRING"
                result.append(" ")
                i += 1
            elif text[i] == "'":
                state = "CHAR"
                result.append(" ")
                i += 1
            else:
                result.append(text[i])
                i += 1
                
        elif state == "LINE_COMMENT":
            if text[i] == '\n':
                state = "NORMAL"
                result.append('\n')
            else:
                result.append(' ')
            i += 1
            
        elif state == "BLOCK_COMMENT":
            if text.startswith("/*", i):
                block_comment_depth += 1
                result.append("  ")
                i += 2
            elif text.startswith("*/", i):
                block_comment_depth -= 1
                result.append("  ")
                i += 2
                if block_comment_depth == 0:
                    state = "NORMAL"
            else:
                if text[i] == '\n':
                    result.append('\n')
                else:
                    result.append(' ')
                i += 1
                
        elif state == "STRING":
            if text[i] == '\\':
                result.append("  ")
                i += 2
            elif text[i] == '"':
                state = "NORMAL"
                result.append(" ")
                i += 1
            else:
                if text[i] == '\n':
                    result.append('\n')
                else:
                    result.append(' ')
                i += 1
                
        elif state == "RAW_STRING":
            end_seq = '"' + ('#' * raw_string_hashes)
            if text.startswith(end_seq, i):
                state = "NORMAL"
                result.append(" " * len(end_seq))
                i += len(end_seq)
            else:
                if text[i] == '\n':
                    result.append('\n')
                else:
                    result.append(' ')
                i += 1
                
        elif state == "CHAR":
            if text[i] == '\\':
                result.append("  ")
                i += 2
            elif text[i] == "'":
                state = "NORMAL"
                result.append(" ")
                i += 1
            else:
                if text[i] == '\n':
                    result.append('\n')
                else:
                    result.append(' ')
                i += 1

    scrubbed_strings = "".join(result)
    
    # Pass 2: Remove attributes (#[...] and #![...])
    result = []
    i = 0
    n = len(scrubbed_strings)
    attr_depth = 0
    
    while i < n:
        if attr_depth == 0:
            if scrubbed_strings.startswith("#[", i):
                attr_depth = 1
                result.append("  ")
                i += 2
            elif scrubbed_strings.startswith("#![", i):
                attr_depth = 1
                result.append("   ")
                i += 3
            else:
                result.append(scrubbed_strings[i])
                i += 1
        else:
            if scrubbed_strings[i] == '[':
                attr_depth += 1
                result.append(" ")
                i += 1
            elif scrubbed_strings[i] == ']':
                attr_depth -= 1
                result.append(" ")
                i += 1
            else:
                if scrubbed_strings[i] == '\n':
                    result.append('\n')
                else:
                    result.append(" ")
                i += 1

    return "".join(result)

class FileAnalyzer:
    """Encapsulates the state for analyzing a single file."""
    
    def __init__(self, path: Path, boundaries: set[str]):
        self.path = path
        self.raw_text = path.read_text(encoding="utf-8")
        self.scrubbed_text = clean_rust_code(self.raw_text)
        self.lines = self.scrubbed_text.splitlines()
        self.raw_lines = self.raw_text.splitlines()
        self.boundaries = boundaries
        self.rel_path = str(path.relative_to(ROOT))
        self.base_module = module_name(path)
        
        self.nested_modules: list[tuple[str, int]] = []
        self.brace_depths = {'{': 0, '(': 0, '[': 0}
        self.macro_body_depth: int | None = None
        self.current_macro_name: str | None = None
        self.macro_open_char: str | None = None
        
        self.in_use_stmt = False
        self.use_stmt_lines = []
        self.use_start_lineno = 0
        self.seen_code = False
        
        self.finding_keys: set[tuple[str, str, int, str]] = set()
        self.edges: list[ImportEdge] = []
        self.alias_imports: list[AliasImport] = []
        self.macro_invocations: list[tuple[str, str, int]] = []
        self.macro_definitions: dict[str, set[str]] = {}
        self.file_imports: set[str] = set()

    def analyze(self) -> None:
        for lineno, line in enumerate(self.lines, start=1):
            self._analyze_line(lineno, line)

    def _analyze_line(self, lineno: int, line: str) -> None:
        is_use_line = False
        if not self.in_use_stmt:
            m = re.match(r"^\s*(?:pub(?:\([^)]+\))?\s+)?use\s+(.*)", line)
            if m:
                self.in_use_stmt = True
                is_use_line = True
                self.use_start_lineno = lineno
                self.use_stmt_lines = [m.group(1)]
                if ";" in line:
                    self.in_use_stmt = False
                    self._process_collected_use_stmt()
        else:
            is_use_line = True
            self.use_stmt_lines.append(line)
            if ";" in line:
                self.in_use_stmt = False
                self._process_collected_use_stmt()

        # Check for modular structure
        mod_match = INLINE_MOD_RE.match(line)
        if mod_match:
            next_depth = self.brace_depths['{'] + line.count("{") - line.count("}")
            self.nested_modules.append((mod_match.group(1), next_depth))

        mr_match = re.search(r"\bmacro_rules!\s+([A-Za-z_][A-Za-z0-9_]*)\s*([{\[(])", line)
        if self.macro_body_depth is None and mr_match:
            self.current_macro_name = mr_match.group(1)
            self.macro_open_char = mr_match.group(2)
            self.macro_body_depth = self.brace_depths[self.macro_open_char]

        inside_macro = (self.macro_body_depth is not None and self.macro_open_char and self.brace_depths[self.macro_open_char] > self.macro_body_depth)
        inside_test_module = any(is_test_module_name(name) for name, _depth in self.nested_modules)
        is_mod_line = line.lstrip().startswith("mod ") or line.lstrip().startswith("pub mod ")

        if is_use_line and self.seen_code and not inside_test_module and not inside_macro and not self.nested_modules:
            self.finding_keys.add(("use_not_at_top", self.rel_path, lineno, "use statement found after other code"))

        if not is_use_line and not inside_test_module and not is_mod_line and not mr_match:
            stripped = line.strip()
            is_ignored_prefix = stripped.startswith("pub(in ") or stripped.startswith("macro_rules!")
            if stripped and stripped != "}" and stripped != "];":
                if not stripped.startswith("extern crate ") and not stripped.startswith("compile_error!") and not is_ignored_prefix:
                    self.seen_code = True
            if not is_ignored_prefix:
                self._handle_code_line(line, lineno, inside_macro)

        self._update_braces(line)

    def _process_collected_use_stmt(self) -> None:
        full_use = " ".join(self.use_stmt_lines).split(";", 1)[0].strip()
        self.use_stmt_lines = []
        
        source_parts = self.base_module.split("::") + [name for name, _depth in self.nested_modules]
        if source_parts == [""]: source_parts = []
        source_mod = "::".join(source_parts)
        source_boundary = boundary_for_module(source_mod, self.boundaries)
        inside_test_module = any(is_test_module_name(name) for name, _depth in self.nested_modules)
        
        if inside_test_module:
            return
            
        self._handle_use_statement(full_use, source_mod, source_boundary, self.use_start_lineno)

    def _handle_use_statement(self, raw_use: str, source_mod: str, source_boundary: str | None, lineno: int) -> None:
        for raw_target, alias in expand_use_tree(raw_use):
            normalized = resolve_relative_path(source_mod, raw_target)
            self.file_imports.add(normalized)
            self.edges.append(
                ImportEdge(
                    source=source_mod,
                    target=normalized,
                    raw=raw_target,
                    path=self.rel_path,
                    line=lineno,
                    internal=is_internal_target(normalized),
                )
            )

            if alias and alias != "_" and not is_whitelisted_alias_import_for_path(self.rel_path, alias, normalized):
                usage_count = count_alias_usage(self.raw_lines, alias, lineno)
                self.alias_imports.append(
                    AliasImport(
                        path=self.rel_path,
                        line=lineno,
                        source_module=source_mod,
                        target=normalized,
                        alias=alias,
                        usage_count=usage_count,
                    )
                )
                if usage_count == 0:
                    self.finding_keys.add(("unused_alias_import", self.rel_path, lineno, f"{alias} -> {normalized}"))

            if raw_target.startswith("super::super"):
                self.finding_keys.add(("deep_relative_import", self.rel_path, lineno, raw_target))

            if raw_target.startswith("crate::"):
                if source_boundary is not None:
                    target_boundary = boundary_for_module(normalized, self.boundaries)
                    if target_boundary == source_boundary:
                        self.finding_keys.add(("self_referential_crate_import", self.rel_path, lineno, f"{raw_target} (should be relative)"))

            if is_shallow_crate_import(raw_target) and not is_whitelisted_crate_root_import(raw_target):
                self.finding_keys.add(("crate_root_import", self.rel_path, lineno, raw_target))

            if "::logic::" in normalized and normalized.startswith("src::"):
                logic_boundary = normalized.split("::logic::", 1)[0]
                if source_boundary != logic_boundary:
                    self.finding_keys.add(("cross_boundary_logic_import", self.rel_path, lineno, f"{source_mod} -> {normalized}"))
                elif not (source_mod == logic_boundary or source_mod.endswith("::api") or source_mod.startswith(f"{logic_boundary}::logic")):
                    after_logic = normalized.split("::logic::", 1)[1]
                    if "::" in after_logic:
                        self.finding_keys.add(("deep_logic_surface_bypass", self.rel_path, lineno, f"{source_mod} -> {normalized}"))
                    else:
                        self.finding_keys.add(("logic_surface_import", self.rel_path, lineno, f"{source_mod} -> {normalized}"))

    def _handle_code_line(self, line: str, lineno: int, inside_macro: bool) -> None:
        for inv_match in MACRO_INVOCATION_RE.findall(line):
            if inv_match not in {"macro_rules"}:
                self.macro_invocations.append((self.rel_path, inv_match, lineno))

        for match in qualified_path_matches(line):
            if inside_macro and self.current_macro_name:
                resolved_path = match
                if resolved_path.startswith("$crate::"):
                    resolved_path = "src::" + resolved_path[8:]
                elif resolved_path.startswith("crate::"):
                    resolved_path = "src::" + resolved_path[7:]
                self.macro_definitions.setdefault(self.current_macro_name, set()).add(resolved_path)
            else:
                self.finding_keys.add(("inline_qualified_path", self.rel_path, lineno, match))

    def _update_braces(self, line: str) -> None:
        self.brace_depths['{'] += line.count("{") - line.count("}")
        self.brace_depths['('] += line.count("(") - line.count(")")
        self.brace_depths['['] += line.count("[") - line.count("]")
        
        if self.macro_body_depth is not None and self.macro_open_char and self.brace_depths[self.macro_open_char] <= self.macro_body_depth:
            self.macro_body_depth = None
            self.current_macro_name = None
            self.macro_open_char = None
            
        while self.nested_modules and self.brace_depths['{'] < self.nested_modules[-1][1]:
            self.nested_modules.pop()


def module_name(path: Path) -> str:
    rel = path.relative_to(ROOT).with_suffix("")
    parts = list(rel.parts)
    if parts and parts[-1] == "mod":
        parts.pop()
    return "::".join(parts)


def staircase_boundaries() -> set[str]:
    boundaries: set[str] = set()
    for directory in SRC.rglob("*"):
        if not directory.is_dir():
            continue
        if (directory / "api.rs").exists() and (directory / "logic").is_dir():
            rel = directory.relative_to(ROOT)
            boundaries.add("::".join(rel.parts))
    return boundaries


def is_test_path(path: Path) -> bool:
    if any(part in TEST_PATH_PARTS for part in path.parts):
        return True
    stem = path.stem
    return stem == "tests" or stem.endswith("_tests")


def is_test_module_name(name: str) -> bool:
    return name == "tests" or name.endswith("_tests")


def split_top_level(text: str, delimiter: str = ",") -> list[str]:
    items: list[str] = []
    current: list[str] = []
    depth = 0
    for char in text:
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
        if char == delimiter and depth == 0:
            item = "".join(current).strip()
            if item:
                items.append(item)
            current = []
            continue
        current.append(char)
    tail = "".join(current).strip()
    if tail:
        items.append(tail)
    return items


def split_alias(path: str) -> str:
    parts = re.split(r"\s+as\s+", path, maxsplit=1)
    return parts[0].strip()


def parse_alias(path: str) -> tuple[str, str | None]:
    parts = re.split(r"\s+as\s+", path, maxsplit=1)
    target = parts[0].strip()
    alias = parts[1].strip() if len(parts) == 2 else None
    return target, alias


def expand_use_tree(tree: str) -> list[tuple[str, str | None]]:
    tree = tree.strip()
    if tree.startswith("{") and tree.endswith("}"):
        paths: list[tuple[str, str | None]] = []
        for item in split_top_level(tree[1:-1]):
            paths.extend(expand_use_tree(item))
        return paths

    depth = 0
    group_start = -1
    prefix_end = -1
    for idx, char in enumerate(tree):
        if char == "{":
            if depth == 0:
                group_start = idx
                prefix_end = idx - 2
            depth += 1
        elif char == "}":
            depth -= 1
            if depth == 0 and group_start != -1:
                prefix = tree[:prefix_end].strip()
                group = tree[group_start + 1 : idx]
                suffix = tree[idx + 1 :].strip()
                if suffix.startswith("::"):
                    suffix = suffix[2:]
                paths: list[tuple[str, str | None]] = []
                for item in split_top_level(group):
                    expanded = expand_use_tree(item)
                    for part, alias in expanded:
                        combined = prefix if part == "self" else f"{prefix}::{part}"
                        if suffix:
                            combined = f"{combined}::{suffix}"
                        paths.append((combined, alias))
                return paths

    target, alias = parse_alias(tree)
    return [(target, alias)]


def resolve_relative_path(source: str, target: str) -> str:
    target = split_alias(target)
    if target == "self":
        return source
    if target == "super":
        return "::".join(source.split("::")[:-1])
    if target.startswith("crate::"):
        return f"src::{target[7:]}"
    if target.startswith("self::"):
        return f"{source}::{target[6:]}"
    if target.startswith("super::"):
        parts = source.split("::")
        rest = target
        while rest.startswith("super::"):
            rest = rest[7:]
            if len(parts) > 1:
                parts.pop()
        if rest:
            parts.extend(rest.split("::"))
        return "::".join(parts)
    return target


def boundary_for_module(module: str, boundaries: set[str]) -> str | None:
    candidates = [
        boundary
        for boundary in boundaries
        if module == boundary or module.startswith(f"{boundary}::")
    ]
    if not candidates:
        return None
    return max(candidates, key=lambda item: item.count("::"))


def is_internal_target(target: str) -> bool:
    return target.startswith("src::")


def is_shallow_crate_import(raw_target: str) -> bool:
    if not raw_target.startswith("crate::"):
        return False
    return "::" not in raw_target[7:]


def is_whitelisted_crate_root_import(raw_target: str) -> bool:
    return raw_target in CONFIG["whitelisted_crate_root_imports"]


def is_whitelisted_alias_import(alias: str, target: str) -> bool:
    return (alias, target) in CONFIG["whitelisted_alias_imports"]


def is_whitelisted_alias_import_for_path(path: str, alias: str, target: str) -> bool:
    if is_whitelisted_alias_import(alias, target):
        return True
    for prefix, allowed in CONFIG["whitelisted_alias_import_prefixes"].items():
        if path.startswith(prefix) and (alias, target) in allowed:
            return True
    return False


def qualified_path_matches(line: str) -> list[str]:
    matches: list[str] = []
    for qp in QUALIFIED_PATH_RE.finditer(line):
        candidate = qp.group(1).strip()
        parts = candidate.split("::")
        
        # SMART CODEBASE ABSTRACTION:
        # Allow paths with exactly one '::' ONLY IF the first segment is PascalCase (a Type) 
        # or a primitive numeric type. This allows idiomatic type-associated paths (e.g. `Scalar::epsilon`) 
        # but bans module-associated paths (e.g. `float_ops::cos`), enforcing `use` for modules.
        if len(parts) == 2:
            head = parts[0]
            if head not in ("crate", "super", "self"):
                is_pascal_case = head and head[0].isupper()
                is_primitive = head in CONFIG["rust_primitive_types"]
                is_whitelisted_module = head in ("float_ops", "int_math", "rational_math", "libm", "alloc")
                if is_pascal_case or is_primitive or is_whitelisted_module:
                    continue
            if candidate in CONFIG.get("whitelisted_inline_qualified_paths", set()):
                continue

        # Skip double-underscore items
        if any(seg.startswith("__") for seg in parts):
            continue
        if re.match(r"^[A-Za-z_][A-Za-z0-9_]*::<", candidate):
            continue
        if "::{" in candidate:
            continue
            
        head = parts[0]
        if head in RUST_KEYWORDS or keyword.iskeyword(head):
            matches.append(candidate)
            continue
            
        matches.append(candidate)
    return matches


def count_alias_usage(lines: list[str], alias: str, start_line: int) -> int:
    pattern = re.compile(rf"\b{re.escape(alias)}\b")
    count = 0
    for line in lines[start_line:]:
        line = re.sub(r"//.*$", "", line)
        count += len(pattern.findall(line))
    return count


def collect_findings() -> tuple[list[Finding], list[ImportEdge], set[str], list[AliasImport]]:
    boundaries = staircase_boundaries()
    finding_keys: set[tuple[str, str, int, str]] = set()
    edges: list[ImportEdge] = []
    alias_imports: list[AliasImport] = []
    
    macro_definitions: dict[str, set[str]] = {}
    macro_invocations: list[tuple[str, str, int]] = []
    file_imports: dict[str, set[str]] = {}

    for path in sorted(SRC.rglob("*.rs")):
        if is_test_path(path):
            continue
            
        analyzer = FileAnalyzer(path, boundaries)
        analyzer.analyze()
        
        finding_keys.update(analyzer.finding_keys)
        edges.extend(analyzer.edges)
        alias_imports.extend(analyzer.alias_imports)
        
        # Merge macro info
        macro_invocations.extend(analyzer.macro_invocations)
        for macro, reqs in analyzer.macro_definitions.items():
            macro_definitions.setdefault(macro, set()).update(reqs)
        
        file_imports[analyzer.rel_path] = analyzer.file_imports

    # Process macros globally
    macro_import_missing_callers: dict[tuple[str, str], list[str]] = {}
    for file_path, inv_macro, lineno in macro_invocations:
        if inv_macro in macro_definitions:
            required_imports = macro_definitions[inv_macro]
            actual_imports = file_imports.get(file_path, set())
            
            for req_import in required_imports:
                caller_mod = file_path.replace('.rs', '').replace('/', '::')
                if caller_mod.endswith('::mod'):
                    caller_mod = caller_mod[:-5]
                req_mod = "::".join(req_import.split("::")[:-1])
                
                is_available = (req_import in actual_imports) or (caller_mod == req_mod)
                if not is_available:
                    macro_import_missing_callers.setdefault((inv_macro, req_import), []).append(f"{file_path}:{lineno}")

    for macro, req_imports in macro_definitions.items():
        for req_import in req_imports:
            missing_callers = macro_import_missing_callers.get((macro, req_import), [])
            if not missing_callers:
                finding_keys.add(
                    ("macro_inline_path_safe_to_remove", f"macro: {macro}", 1,
                     f"Macro uses inline '{req_import}' - SAFE to simplify because all callers import it")
                )

    findings = [Finding(kind=kind, path=path, line=line, text=text) for kind, path, line, text in sorted(finding_keys)]
    alias_imports.sort(key=lambda item: (item.path, item.line, item.alias, item.target))
    return findings, edges, boundaries, alias_imports


def write_dot(edges: list[ImportEdge], boundaries: set[str], dot_path: Path) -> None:
    dot_path.parent.mkdir(parents=True, exist_ok=True)
    nodes = {edge.source for edge in edges}
    nodes.update(edge.target for edge in edges if edge.internal)

    boundary_map: dict[str, list[str]] = {boundary: [] for boundary in sorted(boundaries)}
    loose_nodes: list[str] = []
    
    for node in sorted(nodes):
        boundary = boundary_for_module(node, boundaries)
        if boundary:
            boundary_map.setdefault(boundary, []).append(node)
        else:
            loose_nodes.append(node)

    with dot_path.open("w", encoding="utf-8") as fh:
        fh.write("digraph imports {\n  rankdir=\"LR\";\n  graph [fontname=\"monospace\"];\n  node [shape=\"box\", fontname=\"monospace\"];\n  edge [fontname=\"monospace\"];\n")
        cluster_index = 0
        for boundary, boundary_nodes in boundary_map.items():
            if not boundary_nodes:
                continue
            fh.write(f"  subgraph cluster_{cluster_index} {{\n    label=\"{boundary}\";\n    color=\"lightgrey\";\n")
            for node in sorted(boundary_nodes):
                fill = "white"
                if node == boundary:
                    fill = "lightblue"
                elif "::logic" in node:
                    fill = "lightyellow"
                elif node.endswith("::api"):
                    fill = "honeydew"
                fh.write(f'    "{node}" [style="filled", fillcolor="{fill}"];\n')
            fh.write("  }\n")
            cluster_index += 1

        for node in loose_nodes:
            fh.write(f'  "{node}";\n')

        for edge in edges:
            target = edge.target
            if not edge.internal:
                if "::" in target:
                    target = target.split("::", 1)[0]
                fh.write(f'  "{target}" [shape="ellipse", style="dashed"];\n')
            safe_source = edge.source.replace('"', '\\"')
            safe_target = target.replace('"', '\\"')
            fh.write(f'  "{safe_source}" -> "{safe_target}";\n')
        fh.write("}\n")


def structure_findings(boundaries: set[str]) -> list[Finding]:
    findings: list[Finding] = []
    for boundary in sorted(boundaries):
        rel = ROOT.joinpath(*boundary.split("::"))
        if not (rel / "mod.rs").exists():
            findings.append(Finding(kind="missing_boundary_mod", path=str(rel.relative_to(ROOT)), line=1, text=f"{boundary} has api.rs and logic/ but no mod.rs"))
        if not (rel / "logic" / "mod.rs").exists():
            findings.append(Finding(kind="missing_logic_mod", path=str((rel / "logic").relative_to(ROOT)), line=1, text=f"{boundary}::logic is missing mod.rs"))
    return findings


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--json", action="store_true", help="Emit JSON findings")
    parser.add_argument("--dot", default=str(ROOT / "tools" / "import_graph.dot"), help="Path to write Graphviz DOT output")
    parser.add_argument("--limit", type=int, default=200, help="Maximum number of findings to print in text mode")
    args = parser.parse_args()

    findings, edges, boundaries, alias_imports = collect_findings()
    findings.extend(structure_findings(boundaries))
    findings.sort(key=lambda item: (item.kind, item.path, item.line, item.text))
    write_dot(edges, boundaries, Path(args.dot))
    
    counts = Counter(finding.kind for finding in findings)
    unused_alias_count = sum(1 for item in alias_imports if item.usage_count == 0)

    if args.json:
        print(json.dumps({
            "root": str(ROOT),
            "dot": str(args.dot),
            "boundary_count": len(boundaries),
            "boundaries": sorted(boundaries),
            "counts": dict(counts),
            "alias_import_count": len(alias_imports),
            "unused_alias_import_count": unused_alias_count,
            "alias_imports": [asdict(alias_import) for alias_import in alias_imports],
            "findings": [asdict(finding) for finding in findings],
        }, indent=2))
        return 0

    print(f"DOT graph written to: {args.dot}")
    print(f"Boundaries detected: {len(boundaries)}")
    print(f"Findings: {len(findings)}")
    print(f"Aliased imports: {len(alias_imports)}")
    print(f"Unused aliased imports: {unused_alias_count}")
    for kind, count in sorted(counts.items()):
        print(f"  {kind}: {count}")
    for finding in findings[: args.limit]:
        print(f"{finding.kind}: {finding.path}:{finding.line}: {finding.text}")
        
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
