use anyhow::Result;
use serde::Serialize;
use std::fs;
use std::io::Write;
use std::path::Path;

pub fn ensure_parent_dir(path: &Path) -> Result<()> {
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)?;
    }
    Ok(())
}

pub fn copy_file(source: &Path, destination: &Path) -> Result<u64> {
    ensure_parent_dir(destination)?;
    Ok(fs::copy(source, destination)?)
}

pub fn link_or_copy_file(source: &Path, destination: &Path) -> Result<()> {
    ensure_parent_dir(destination)?;
    if destination.exists() {
        fs::remove_file(destination)?;
    }
    match fs::hard_link(source, destination) {
        Ok(()) => Ok(()),
        Err(_) => {
            fs::copy(source, destination)?;
            Ok(())
        }
    }
}

pub fn write_json<T: Serialize>(path: &Path, value: &T) -> Result<()> {
    ensure_parent_dir(path)?;
    let file = fs::File::create(path)?;
    serde_json::to_writer_pretty(file, value)?;
    Ok(())
}

pub fn write_toml<T: Serialize>(path: &Path, value: &T) -> Result<()> {
    ensure_parent_dir(path)?;
    let mut file = fs::File::create(path)?;
    let text = toml::to_string_pretty(value)?;
    file.write_all(text.as_bytes())?;
    file.flush()?;
    Ok(())
}
