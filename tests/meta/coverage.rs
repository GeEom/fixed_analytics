//! Ensures all src files have corresponding test files

#[cfg(test)]
mod tests {
    use std::collections::HashSet;
    use std::fs;
    use std::io;
    use std::path::Path;

    #[test]
    fn test_all_src_files_have_unit_tests() {
        let src_dir = Path::new("src");
        let tests_dir = Path::new("tests/unit");

        let src_paths = collect_relative_paths(src_dir, src_dir).unwrap_or_default();
        let test_paths = if tests_dir.exists() {
            collect_relative_paths(tests_dir, tests_dir).unwrap_or_default()
        } else {
            HashSet::new()
        };

        let mut missing_tests = Vec::new();

        for src_path in &src_paths {
            // Skip entry points and module files
            if src_path == "main.rs" || src_path == "lib.rs" || src_path.ends_with("mod.rs") {
                continue;
            }

            if !test_paths.contains(src_path) {
                missing_tests.push(src_path.clone());
            }
        }

        assert!(
            missing_tests.is_empty(),
            "Missing unit test files:\n{}",
            missing_tests
                .iter()
                .map(|p| format!("  - src/{p} -> tests/unit/{p}"))
                .collect::<Vec<_>>()
                .join("\n")
        );
    }

    fn collect_relative_paths(dir: &Path, base: &Path) -> Result<HashSet<String>, io::Error> {
        let mut paths = HashSet::new();
        if dir.is_dir() {
            for entry in fs::read_dir(dir)? {
                let path = entry?.path();
                let relative = path
                    .strip_prefix(base)
                    .map_err(|_| io::Error::other("strip prefix failed"))?
                    .to_string_lossy()
                    .to_string();

                if path.is_dir() {
                    paths.insert(relative.clone());
                    paths.extend(collect_relative_paths(&path, base)?);
                } else if path.extension().and_then(|e| e.to_str()) == Some("rs") {
                    paths.insert(relative);
                }
            }
        }
        Ok(paths)
    }
}
