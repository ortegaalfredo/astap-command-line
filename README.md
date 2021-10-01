# astap-cli

ASTAP (Astrometric Stacker Program) command line version, with commands specified through a JSON file in order to make it easily integrable in other applications.

Example of stacked image:

![Output.jpg](https://github.com/ortegaalfredo/astap-command-line/blob/main/Samples/Output.jpg?raw=true)


This is a branch of ASTAP (https://www.hnsky.org/astap.htm) command-line utility, with additional support for sigma-clip stacking and all needed processing.

ASTAP already provides a command-line version but only capable of solving. This project adds creation of master darks, master flats and sigma-clipping stacking with astrometric or star-based alignment.

## Building

This project requires Lazarus 2.2.0 and FreePascal 3.2.2.
Binaries for win32 and static linux (Built on ubuntu 20.04.3 LTS) are provided.

## Calibration

It's optional but recommendable that you calibrate the images before stacking and solving. This must be done using the traditional Dark/Flat and DarkFlats images. You can use this utility to average all calibration images in the following way:

### 1) Generate master darks

The command to generate master darks is:

```
> astap-cli -stack 1-create-master-dark.json
```

And in 1-create-master-dark.json all darks images must be specified in this way:

```
{
"Command" : "CreateMasterDark",
"Darks": [
	  "Darks\\bin1-120s_2188.fit",
	  "Darks\\bin1-120s_2189.fit"
          ],
"Output": "Darks\\MasterDark.fit"
}
```

### 2) Generate master flats

The command to generate master flats is:

```
> astap-cli -stack 1-create-master-darks.json
```

And in 1-create-master-flat.json all darks images must be specified in this way:

```
{
"Command" : "CreateMasterFlat",
"DarkFlats": [
	  "Darks\\bin1-120s_2188.fit",
	  "Darks\\bin1-120s_2189.fit"
          ],
"Flats": [
	  "Flats\\bin1_1983.fit",
	  "Flats\\bin1_1984.fit"
          ],
"Output": "Flats\\MasterFlat.fit"
}
```

Note that this command differs from the creation of the Master Dark in that you can optionally specify DarkFlats images that will be averaged and used to correct the Master Flat.

## Sigma-clipping Stacking

ASTAP provides several stacking algorithms but currently astap-cli.exe only supports sigma-clipping. Also, this utility provides two alignment methods:
1) Star Alignment: Compare stars in pictures to align them.
2) Astrometric Alignment: Search the stars in a star database to align the images.

The command to stack a series of images is:
```
> astap-cli -stack C:\astro\astap\command-line_version\stacktest.json

```

This is the format of the JSON file:

```
{
"Command" : "Stack",
"Align"   : "stars",
"MasterDark" :     "Darks\\MasterDark.fit",
"MasterFlat" :     "Flats\\MasterFlat.fit",
"Lights": [
	  "Lights\\helix-180-HA_4628.fit",
	  "Lights\\helix-180-HA_4629.fit"
          ],
"Output": "Lights\\Output.fit"
}
```

The "Align" parameter can be either "stars" or "astrometry" to select each alignment method.
The "MasterDark" and "MasterFlat" specify the result image of processing Dark and Flat images.
The "Lights" parameter specify the list of images to stack.
The "Output" parameter specify the output file.

Note that in the same way as ASTAP, you can specify command line options to modify parameters. For example, to specify the position of the star-database when using astrometry solving, the command line might look like this:

```
> astap-cli -fov 1 -d "C:\Program Files\astap" -stack 3-stack.json
```

Where -fov specify the approximate FOV of the image in degrees (This parameter is essential to achieve a successful astrometric solution) and -d parameter specify the star-database installation path. Currently only the H18 star database is supported (Download from https://sourceforge.net/projects/astap-program/files/star_databases/h18_star_database_mag18_astap.exe/download)

## Other command-line options

Command line options are the same as in the regular ASTAP command-line tool. Options affects the astrometric solver, but not the stacking, as it's specified completely using JSON files.

| command  |      Parameter      |  unit | Remarks
|----------|:-------------------:|------:|------------:|
| -h | | | help info
| -help |       |   | help info
| -f | file_name |     | File to solve astrometric.
| -f | stdin     |     | Will accept raw images via stdin. Raw format according INDI.
|-r | radius_search_field |	degrees|The program will search in a square spiral around the start position up to this radius *
|-fov| 	field_height_of_image |	degrees |	Optional. Normally calculated from FITS header. Use value 0 for auto.  If 0 is specified the fov found by solving it will be saved for next time.(learn mode)  *
|-ra |	center_right_ascension | hours | Optional start value. Normally calculated from FITS header.
|-spd |	center_south_pole_distance (dec+90) |	degrees |	Normally calculated from FITS header * The declination is given in south pole distance, so always positive.
|-z 	| down_sample_factor |0,1,2,3,4 | 	Down sample prior to solving. Also called binning. A value "0" will result in auto selection downsampling. *
|-s |	max_number_of_stars | |	Limits the number of star used for the solution. Typical value 500. *
|-t |	tolerance |  | Tolerance used to compare quads. Typical value  0.007. *
|-m | 	minimum star size |	arcsec	| This could be used to filter out hot pixels.
|-speed |	mode | slow / auto | 	"slow" is forcing of reading a larger area from the star database (more overlap)  to improve detection. *
|-o | file |	Name the output files with this base path & file name
|-d | path |	Specify a path to the star database
|-analyse  |snr_minimum | |	Analyse only and report HFD. Windows: errorlevel is the median HFD * 100M + number of stars used. So the HFD is trunc(errorlevel/1M)/100. For Linux and macOS the info is send to stdout only.
|-analyse2 | snr_minimum |		|Both analyse and solve the image
| -extract | snr_minimum |		| As analyse option but additionally write a .csv file with the detected stars info.
|-sqm	|pedestal |	|	measure the sky background value in magn/arcsec2    relative to the stars. The pedestal is the mean value of a dark.
|-focus1 |	file1.fits -focus2 file2.fits .......	|	| Find best focus point for four or more images using curve fitting. Windows: errorlevel is focuspos*1E4 + rem.error*1E3. Linux: see stdout
|-annotate	|	|	| Produce a deep sky annotated jpeg file with same name as input file extended with _annotated *
|-debug		|	|   | Show GUI and stop prior to solving
|-log		|	|   |Write solver log to file with extension .log
| -tofits	|binning|	1,2,3,4	| Produce binned FITS file from input png/jpg *
| -update   |	|   | Update the fits header with the found solution *
|-wcs       |   |   |	Write a .wcs file  in similar format as Astrometry.net. Else text style.


