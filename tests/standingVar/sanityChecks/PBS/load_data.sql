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
.import /g/data/ht96/nb9894/standingVar/sanityChecks/slim_muts.csv slim_muts

ALTER TABLE slim_muts
DROP COLUMN const
DROP COLUMN chi;