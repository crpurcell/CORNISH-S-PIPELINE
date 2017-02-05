
CREATE TABLE observations (
timeStamp VARCHAR(20) PRIMARY KEY, 
dayName VARCHAR(20), 
pntName a VARCHAR(10),
calCode VARCHAR(3),
nAnts INTEGER,
fieldName VARCHAR(50),
RA_deg DOUBLE,
Dec_deg DOUBLE,
goodvis DOUBLE,
badvis DOUBLE,
splitFlag INTEGER,
useFlag INTEGER
);

CREATE TABLE mosaic_obs (
dayPntID VARCHAR(30) PRIMARY KEY,
dayName VARCHAR(20),
pntName VARCHAR(10),
calCode VARCHAR(3),
RA_deg DOUBLE,
Dec_deg DOUBLE,
splitFlag INTEGER,
useFlag INTEGER
);

CREATE TABLE field_coords (
fieldName VARCHAR(50) PRIMARY KEY, 
RA_deg DOUBLE,
Dec_deg DOUBLE,
calCode VARCHAR(3));

CREATE TABLE pointings (
pntName VARCHAR(10),
calCode VARCHAR(3),
RA_deg DOUBLE,
Dec_deg DOUBLE,
fieldName VARCHAR(50));

CREATE TABLE field_images (
fieldName VARCHAR(50) PRIMARY KEY,
tileID INTEGER,
RA_deg DOUBLE,
Dec_deg DOUBLE,
no_fields INTEGER,
pix_max DOUBLE,
pix_min DOUBLE,
npix_x INTEGER,
npix_y INTEGER,
pixscale_x DOUBLE,
pixscale_y DOUBLE,
bm_maj DOUBLE,
bm_min DOUBLE,
bm_pa DOUBLE,
rms_quality DOUBLE,
rms_mult DOUBLE,
uvrange_min DOUBLE,
uvrange_max DOUBLE,
obs_includedVARCHAR(50));
