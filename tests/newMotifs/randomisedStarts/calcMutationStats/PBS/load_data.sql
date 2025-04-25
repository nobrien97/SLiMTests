PRAGMA journal_mode = OFF;
PRAGMA synchronous = OFF;
PRAGMA cache_size = 7500000;

CREATE TABLE tab_epistasis(
    gen             INTEGER,
    seed            INTEGER,
    modelindex      INTEGER,
    mutType_ab      TEXT,
    wa              REAL,
    wb              REAL,
    wab             REAL,
    wwt             REAL,
    ew              REAL,
    ew_s            REAL
);

.mode csv
.import /g/data/ht96/nb9894/newMotifs/randomisedStarts/calcMutationStats/d_epistasis.csv tab_epistasis

CREATE TABLE tab_epistasis_freq(
    gen             INTEGER,
    seed            INTEGER,
    modelindex      INTEGER,
    mutType_ab      TEXT,
    wa              REAL,
    wb              REAL,
    wab             REAL,
    wwt             REAL,
    ew              REAL,
    ew_s            REAL
);

.mode csv
.import /g/data/ht96/nb9894/newMotifs/calcMutationStats/d_epistasis_freqweight.csv tab_epistasis_freq
