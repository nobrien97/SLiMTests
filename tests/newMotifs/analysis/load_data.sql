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
    originGen     INTEGER,
    effect        REAL,
    freq          REAL,
    count         INTEGER,
    fixGen        INTEGER
);

.mode csv
.import /g/data/ht96/nb9894/newMotifs/slim_muts.csv slim_muts
.import /g/data/ht96/nb9894/newMotifs/extraReplicates/slim_muts.csv slim_muts

ALTER TABLE slim_muts
DROP COLUMN const
DROP COLUMN chi;

CREATE TABLE slim_qg(
    gen             INTEGER,
    seed            INTEGER,
    modelindex      INTEGER,
    meanH           REAL,
    trait1_mean     REAL,
    trait2_mean     REAL,
    trait3_mean     REAL,
    trait4_mean     REAL,
    trait1_var      REAL,
    trait2_var      REAL,
    trait3_var      REAL,
    trait4_var      REAL,
    dist            REAL,
    dist1           REAL,
    dist2           REAL,
    dist3           REAL,
    dist4           REAL,
    mean_w          REAL,
    var_w           REAL,
    deltaPheno      REAL,
    deltaW          REAL,
    meanMC1         REAL,
    meanMC2         REAL,
    meanMC3         REAL,
    meanMC4         REAL,
    meanMC5         REAL,
    meanMC6         REAL,
    meanMC7         REAL,
    meanMC8         REAL,
    meanMC9         REAL,
    meanMC10        REAL,
    meanMC11        REAL
);

.mode csv
.import /g/data/ht96/nb9894/newMotifs/slim_qg.csv slim_qg

.backup /iointensive/standingVarMuts.db.bak