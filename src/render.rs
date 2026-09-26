//! SVG -> PNG with resvg, using only the embedded DejaVu Sans (matplotlib's
//! default face), so output never depends on the host's fonts.

use std::{
    fs,
    path::{Path, PathBuf},
    sync::{Arc, OnceLock},
};

use anyhow::{Context, Result};
use resvg::{tiny_skia, usvg};

use crate::Format;

const REGULAR: &[u8] = include_bytes!("../assets/fonts/DejaVuSans.ttf");
const BOLD: &[u8] = include_bytes!("../assets/fonts/DejaVuSans-Bold.ttf");
pub const FONT: &str = "DejaVu Sans";

/// PNG pixels per SVG unit: figures are laid out at 1x and rasterised at 2x.
const SCALE: f32 = 2.0;

fn fonts() -> Arc<usvg::fontdb::Database> {
    static DB: OnceLock<Arc<usvg::fontdb::Database>> = OnceLock::new();
    DB.get_or_init(|| {
        let mut db = usvg::fontdb::Database::new();
        db.load_font_data(REGULAR.to_vec());
        db.load_font_data(BOLD.to_vec());
        db.set_sans_serif_family(FONT);
        Arc::new(db)
    })
    .clone()
}

pub fn png(svg: &str) -> Result<Vec<u8>> {
    let opts = usvg::Options {
        fontdb: fonts(),
        font_family: FONT.into(),
        ..Default::default()
    };
    let tree = usvg::Tree::from_str(svg, &opts).context("parsing generated SVG")?;
    let size = tree
        .size()
        .to_int_size()
        .scale_by(SCALE)
        .context("empty figure")?;
    let mut pixmap =
        tiny_skia::Pixmap::new(size.width(), size.height()).context("figure too large")?;
    pixmap.fill(tiny_skia::Color::WHITE);
    resvg::render(
        &tree,
        tiny_skia::Transform::from_scale(SCALE, SCALE),
        &mut pixmap.as_mut(),
    );
    pixmap.encode_png().context("encoding PNG")
}

/// Write `<stem>.png` and/or `<stem>.svg`; returns the paths written.
pub fn write(svg: &str, stem: &Path, format: Format) -> Result<Vec<PathBuf>> {
    if let Some(dir) = stem.parent() {
        fs::create_dir_all(dir).with_context(|| format!("creating {}", dir.display()))?;
    }
    let with = |ext: &str| {
        let mut s = stem.as_os_str().to_owned();
        s.push(ext);
        PathBuf::from(s)
    };
    let mut out = Vec::new();
    if matches!(format, Format::Png | Format::Both) {
        let p = with(".png");
        fs::write(&p, png(svg)?).with_context(|| format!("writing {}", p.display()))?;
        out.push(p);
    }
    if matches!(format, Format::Svg | Format::Both) {
        let p = with(".svg");
        fs::write(&p, svg).with_context(|| format!("writing {}", p.display()))?;
        out.push(p);
    }
    Ok(out)
}
