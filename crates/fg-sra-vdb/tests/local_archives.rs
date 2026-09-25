//! Tests of path handling that need no network.

use fg_sra_vdb::error::VdbError;
use fg_sra_vdb::manager::VdbManager;

#[test]
fn opening_a_path_containing_a_percent_sign_is_refused() {
    let dir = std::env::temp_dir().join(format!("fg_sra_vdb_{}_percent", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("run%s%s%s%d.sra");
    let mgr = VdbManager::make_read().unwrap();
    assert_eq!(mgr.open_db_read(path.to_str().unwrap()).err(), Some(VdbError::PercentInPath));
}
