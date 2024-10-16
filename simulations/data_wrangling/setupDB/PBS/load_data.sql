PRAGMA journal_mode = OFF;
PRAGMA synchronous = OFF;
PRAGMA cache_size = 7500000;

CREATE TABLE slim_muts(
    gen           INTEGER,
    seed          INTEGER,
    modelindex    INTEGER,
    mutType       INTEGER,
    mutID         INTEGER,
    pos           INTEGER,
    const         INTEGER,
    originGen     INTEGER,
    effect        REAL,
    chi           REAL,
    freq          REAL,
    count         INTEGER,
    fixGen        INTEGER
);

.mode csv
.import /storage_path/runSims/slim_muts.csv slim_muts

ALTER TABLE slim_muts
DROP COLUMN const
DROP COLUMN chi;

CREATE TABLE slim_qg(
    gen             INTEGER,
    seed            INTEGER,
    modelindex      INTEGER,
    meanH           REAL,
    VA              REAL,
    phenomean       REAL,
    phenovar        REAL,
    dist            REAL,
    w               REAL,
    deltaPheno      REAL,
    deltaw          REAL,
    aZ              REAL,
    bZ              REAL,
    KZ              REAL,
    KXZ             REAL       
);

.mode csv
.import /storage_path/runSims/slim_qg.csv slim_qg

.backup /storage_path/runSims/standingVarMuts.db.bak