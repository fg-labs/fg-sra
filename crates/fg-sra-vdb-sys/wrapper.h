/**
 * Wrapper header for bindgen FFI generation.
 *
 * Includes only the VDB C API functions needed by fg-sra-vdb.
 * bindgen processes this single file to produce Rust FFI bindings.
 */

/* Core types and return codes */
#include <klib/defs.h>
#include <klib/rc.h>
#include <klib/namelist.h>

/* VDB Manager, Database, Table, Cursor */
#include <vdb/manager.h>
#include <vdb/database.h>
#include <vdb/table.h>
#include <vdb/cursor.h>
#include <vdb/vdb-priv.h>

/* KDB Metadata */
#include <kdb/meta.h>
#include <kdb/namelist.h>

/* Alignment: references, iterators */
#include <align/manager.h>
#include <align/reference.h>
#include <align/iterator.h>

/* INSDC coordinate types */
#include <insdc/insdc.h>
