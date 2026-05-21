#[cfg(test)]
use anyhow::{Context, Result, bail};
#[cfg(test)]
use std::fs;
#[cfg(test)]
use std::path::Path;
#[cfg(test)]
use std::process::Command;

#[cfg(test)]
pub fn extract_zip_entry(archive: &Path, entry: &str, destination: &Path) -> Result<()> {
    if let Some(parent) = destination.parent() {
        fs::create_dir_all(parent)?;
    }

    let output = Command::new("unzip")
        .arg("-p")
        .arg(archive)
        .arg(entry)
        .output()
        .with_context(|| {
            format!(
                "failed to extract {} from {} using unzip",
                entry,
                archive.display()
            )
        })?;

    if !output.status.success() {
        bail!(
            "unzip failed for {} from {}: {}",
            entry,
            archive.display(),
            String::from_utf8_lossy(&output.stderr)
        );
    }

    fs::write(destination, output.stdout)?;
    Ok(())
}
