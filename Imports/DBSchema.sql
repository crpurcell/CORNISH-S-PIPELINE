
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
chan INT,
lineStr VARCHAR(50),
stokes VARCHAR(2),
rms DOUBLE,
rmsV DOUBLE,
cutoff DOUBLE,
imsize_pix INT,
cellsize_asec DOUBLE,
robust DOUBLE,
niter INT,
minpatch INT,
speed DOUBLE,
imgClnData VARCHAR(50),
beamMajInt_asec DOUBLE,
beamMinInt_asec DOUBLE,
beamPAint_asec DOUBLE,
beamMaj_asec DOUBLE,
beamMin_asec DOUBLE,
beamPA_asec DOUBLE,
success INT);
