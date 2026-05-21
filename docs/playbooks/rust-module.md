# Rust Module Development Playbook

## Purpose
Standardized way to create new Rust modules in the LineaBact project.

## Rules

1. **Module Structure**
   - Create module in appropriate directory (`src/linear/`, `src/graph/`, etc.)
   - Export public items in `mod.rs`
   - Keep modules focused and cohesive

2. **Error Handling**
   - Use `thiserror` to define custom error types
   - Use `anyhow::Result` in binary/entry functions
   - Always propagate errors instead of using `unwrap()` in production code

3. **Documentation**
   - Add `//!` module-level documentation
   - Document all public functions and important private ones
   - Explain the purpose of the module

4. **Testing**
   - Create corresponding tests in the same file or `tests/` directory
   - For complex logic, add unit tests + at least one integration-style test

5. **Linear Module Specific Rules**
   - When working in `src/linear/`, strictly separate logic for left and right ends
   - Use clear naming: `left_*` and `right_*` functions when appropriate
   - Terminal decisions must be traceable

6. **General Style**
   - Follow Rust 2024 edition idioms
   - Prefer clarity over premature optimization
   - Add TODO comments for future improvements instead of leaving unclear code
