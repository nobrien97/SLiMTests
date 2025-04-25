PRAGMA journal_mode = OFF;
PRAGMA synchronous = OFF;
PRAGMA cache_size = 7500000;

CREATE TABLE tab_epistasis_new_epi(
    gen             INTEGER,
    seed            INTEGER,
    modelindex      INTEGER,
    mutType_ab      TEXT,
    wa              REAL,
    wb              REAL,
    wab             REAL,
    Pwt             REAL,
    Pa              REAL,
    Pb              REAL,
    Pab             REAL,
    ew              REAL,
    ep              REAL
);

.mode csv
.import /g/data/ht96/nb9894/standingVar/calcMutationStats/d_epistasis_s.csv tab_epistasis_new_epi

CREATE TABLE tab_epistasis_freq_new_epi(
    gen             INTEGER,
    seed            INTEGER,
    modelindex      INTEGER,
    mutType_ab      TEXT,
    wa              REAL,
    wb              REAL,
    wab             REAL,
    Pwt             REAL,
    Pa              REAL,
    Pb              REAL,
    Pab             REAL,
    ew              REAL,
    ep              REAL
);

.mode csv
.import /g/data/ht96/nb9894/standingVar/calcMutationStats/d_epistasis_s_freqweight.csv tab_epistasis_freq_new_epi
