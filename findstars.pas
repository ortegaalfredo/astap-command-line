unit FindStars;

{$mode ObjFPC}{$H+}

interface

uses
  Classes, SysUtils,unit_command_line_general,unit_star_align;

  {bin, measure background, find stars}
  procedure bin_and_find_stars(img :image_array;binning:integer;cropping,hfd_min:double;get_hist{update hist}:boolean; out starlist3:star_list; out short_warning : string);

implementation



function read_stars(telescope_ra,telescope_dec,search_field : double; nrstars_required: integer; out nrstars:integer): boolean;{read star from star database}
var
   Bp_Rp, ra2,dec2,
   frac1,frac2,frac3,frac4  : double;
   area1,area2,area3,area4,nrstars_required2  : integer;
begin
  result:=false;{assume failure}
  nrstars:=0;{set counters at zero}
  ra2:=0; {define ra2 value. Prevent ra2 = -nan(0xffffffffffde9) run time failure when first header record is read}

  SetLength(starlist1,2,nrstars_required);{set array length}

  {Assume the search field is at a crossing of four tiles. The search field area, by definition 100% is split in 8%, 15%, 20%, 57% area for each tile.
   There are 500 stars required. It will then retrieve 8% x 500, 15% x 500, 20% x 500, 57% x 500 stars from each tile under the condition these stars are within the green area.
   This will work assuming the star density within the green area is reasonable homogene.}
  find_areas( telescope_ra,telescope_dec, search_field,{var} area1,area2,area3,area4, frac1,frac2,frac3,frac4);{find up to four star database areas for the square image}

  {read 1th area}
  if area1<>0 then {read 1th area}
  begin
    if open_database(telescope_dec,area1)=false then
      exit;{open database file or reset buffer}
    nrstars_required2:=min(nrstars_required,trunc(nrstars_required * frac1));
    while ((nrstars<nrstars_required2) and (readdatabase290(telescope_ra,telescope_dec, search_field, {var} ra2,dec2, mag2,Bp_Rp)) ) do{star 290 file database read. Read up to nrstars_required}
    begin {add star}
      equatorial_standard(telescope_ra,telescope_dec,ra2,dec2,1,starlist1[0,nrstars]{x},starlist1[1,nrstars]{y});{store star CCD x,y position}
      inc(nrstars);
    end;
  end;

  if area2<>0 then {read 2th area}
  begin
    if open_database(telescope_dec,area2)=false then
      exit; {open database file or reset buffer}
    nrstars_required2:=min(nrstars_required,trunc(nrstars_required * (frac1+frac2)));{prevent round up errors resulting in error starlist1}
    while ((nrstars<nrstars_required2) and (readdatabase290(telescope_ra,telescope_dec, search_field, {var} ra2,dec2, mag2,Bp_Rp)) ) do{star 290 file database read. Read up to nrstars_required}
    begin {add star}
      equatorial_standard(telescope_ra,telescope_dec,ra2,dec2,1,starlist1[0,nrstars]{x},starlist1[1,nrstars]{y});{store star CCD x,y position}
      inc(nrstars);
    end;
  end;

  if area3<>0 then {read 3th area}
  begin
    if open_database(telescope_dec,area3)=false then
      exit; {open database file or reset buffer}
    nrstars_required2:=min(nrstars_required,trunc(nrstars_required * (frac1+frac2+frac3)));
    while ((nrstars<nrstars_required2) and (readdatabase290(telescope_ra,telescope_dec, search_field, {var} ra2,dec2, mag2,Bp_Rp)) ) do{star 290 file database read. Read up to nrstars_required}
    begin {add star}
      equatorial_standard(telescope_ra,telescope_dec,ra2,dec2,1,starlist1[0,nrstars]{x},starlist1[1,nrstars]{y});{store star CCD x,y position}
      inc(nrstars);
    end;
  end;

  if area4<>0 then {read 4th area}
  begin
    if open_database(telescope_dec,area4)=false then
     exit; {open database file}
    nrstars_required2:=min(nrstars_required,trunc(nrstars_required * (frac1+frac2+frac3+frac4)));
    while ((nrstars<nrstars_required2) and (readdatabase290(telescope_ra,telescope_dec, search_field, {var} ra2,dec2, mag2,Bp_Rp)) ) do{star 290 file database read. Read up to nrstars_required}
    begin {add star}
      equatorial_standard(telescope_ra,telescope_dec,ra2,dec2,1,starlist1[0,nrstars]{x},starlist1[1,nrstars]{y});{store star CCD x,y position}
      inc(nrstars);
    end;
  end;

//  memo2_message('testareas'+#9+floattostr4(telescope_ra*12/pi)+#9+floattostr4(telescope_dec*180/pi)+#9+inttostr(maga)+#9+inttostr(magb)+#9+inttostr(magc)+#9+inttostr(magd)+#9+floattostr4(frac1)+#9+floattostr4(frac2)+#9+floattostr4(frac3)+#9+floattostr4(frac4)+#9+inttostr(area1)+#9+inttostr(area2)+#9+inttostr(area3)+#9+inttostr(area4));

  if nrstars<nrstars_required then
       SetLength(starlist1,2,nrstars); {fix array length on data for case less stars are found}
  result:=true;{no errors}
end;


procedure binX1_crop(crop {0..1}:double; img : image_array; var img2: image_array);{crop image, make mono, no binning}
  var fitsX,fitsY,k, w,h,  shiftX,shiftY: integer;
      val       : single;
begin
  w:=trunc(crop*width2);  {cropped}
  h:=trunc(crop*height2);

  setlength(img2,1,w,h); {set length of image array}

  shiftX:=round(width2*(1-crop)/2); {crop is 0.9, shift is 0.05*width2}
  shiftY:=round(height2*(1-crop)/2); {crop is 0.9, start at 0.05*height2}

  for fitsY:=0 to h-1 do
    for fitsX:=0 to w-1  do
    begin
      val:=0;
      for k:=0 to naxis3-1 do {all colors and make mono}
         val:=val + img[k,shiftX+fitsx   ,shiftY+fitsY];
      img2[0,fitsX,fitsY]:=val/naxis3;
    end;
  width2:=w;
  height2:=h;
  naxis3:=1;
end;


procedure binX2_crop(crop {0..1}:double; img : image_array; var img2: image_array);{combine values of 4 pixels and crop is required, Result is mono}
  var fitsX,fitsY,k, w,h,  shiftX,shiftY,nrcolors,width5,height5: integer;
      val       : single;
begin
   nrcolors:=Length(img);
   width5:=Length(img[0]);    {width}
   height5:=Length(img[0][0]); {height}

   w:=trunc(crop*width5/2);  {half size & cropped. Use trunc for image 1391 pixels wide like M27 test image. Otherwise exception error}
   h:=trunc(crop*height5/2);

   setlength(img2,1,w,h); {set length of image array}

   shiftX:=round(width5*(1-crop)/2); {crop is 0.9, shift is 0.05*width2}
   shiftY:=round(height5*(1-crop)/2); {crop is 0.9, start at 0.05*height2}

   for fitsY:=0 to h-1 do
      for fitsX:=0 to w-1  do
     begin
       val:=0;
       for k:=0 to nrcolors-1 do {all colors}
         val:=val+(img[k,shiftX+fitsx*2   ,shiftY+fitsY*2]+
                   img[k,shiftX+fitsx*2 +1,shiftY+fitsY*2]+
                   img[k,shiftX+fitsx*2   ,shiftY+fitsY*2+1]+
                   img[k,shiftX+fitsx*2 +1,shiftY+fitsY*2+1])/4;
       img2[0,fitsX,fitsY]:=val/nrcolors;
     end;

   width2:=w;
   height2:=h;
   naxis3:=1;
 end;

procedure binX3_crop(crop {0..1}:double; img : image_array; var img2: image_array);{combine values of 9 pixels and crop is required. Result is mono}
  var fitsX,fitsY,k, w,h,  shiftX,shiftY,nrcolors,width5,height5: integer;
      val       : single;
begin
  nrcolors:=Length(img);
  width5:=Length(img[0]);    {width}
  height5:=Length(img[0][0]); {height}

  w:=trunc(crop*width5/3);  {1/3 size and cropped}
  h:=trunc(crop*height5/3);

  setlength(img2,1,w,h); {set length of image array}

  shiftX:=round(width5*(1-crop)/2); {crop is 0.9, shift is 0.05*width2}
  shiftY:=round(height5*(1-crop)/2); {crop is 0.9, start at 0.05*height2}

  for fitsY:=0 to h-1 do {bin & mono image}
    for fitsX:=0 to w-1  do
    begin
      val:=0;
      for k:=0 to nrcolors-1 do {all colors}
                     val:=val+(img[k,shiftX+fitsX*3   ,shiftY+fitsY*3  ]+
                               img[k,shiftX+fitsX*3   ,shiftY+fitsY*3+1]+
                               img[k,shiftX+fitsX*3   ,shiftY+fitsY*3+2]+
                               img[k,shiftX+fitsX*3 +1,shiftY+fitsY*3  ]+
                               img[k,shiftX+fitsX*3 +1,shiftY+fitsY*3+1]+
                               img[k,shiftX+fitsX*3 +1,shiftY+fitsY*3+2]+
                               img[k,shiftX+fitsX*3 +2,shiftY+fitsY*3  ]+
                               img[k,shiftX+fitsX*3 +2,shiftY+fitsY*3+1]+
                               img[k,shiftX+fitsX*3 +2,shiftY+fitsY*3+2])/9;
       img2[0,fitsX,fitsY]:=val/nrcolors;
    end;
  width2:=w;
  height2:=h;
  naxis3:=1;
end;


procedure binX4_crop(crop {0..1}:double;img : image_array; var img2: image_array);{combine values of 16 pixels and crop is required. Result is mono}
  var fitsX,fitsY,k, w,h,  shiftX,shiftY,nrcolors,width5,height5: integer;
      val       : single;
begin
  nrcolors:=Length(img);
  width5:=Length(img[0]);    {width}
  height5:=Length(img[0][0]); {height}

  w:=trunc(crop*width5/4);  {1/4 size and cropped}
  h:=trunc(crop*height5/4);

  setlength(img2,1,w,h); {set length of image array}

  shiftX:=round(width5*(1-crop)/2); {crop is 0.9, shift is 0.05*width2}
  shiftY:=round(height5*(1-crop)/2); {crop is 0.9, start at 0.05*height2}

  for fitsY:=0 to h-1 do {bin & mono image}
    for fitsX:=0 to w-1  do
    begin
      val:=0;
      for k:=0 to nrcolors-1 do {all colors}
                     val:=val+(img[k,shiftX+fitsX*4   ,shiftY+fitsY*4  ]+
                               img[k,shiftX+fitsX*4   ,shiftY+fitsY*4+1]+
                               img[k,shiftX+fitsX*4   ,shiftY+fitsY*4+2]+
                               img[k,shiftX+fitsX*4   ,shiftY+fitsY*4+3]+
                               img[k,shiftX+fitsX*4 +1,shiftY+fitsY*4  ]+
                               img[k,shiftX+fitsX*4 +1,shiftY+fitsY*4+1]+
                               img[k,shiftX+fitsX*4 +1,shiftY+fitsY*4+2]+
                               img[k,shiftX+fitsX*4 +1,shiftY+fitsY*4+3]+
                               img[k,shiftX+fitsX*4 +2,shiftY+fitsY*4  ]+
                               img[k,shiftX+fitsX*4 +2,shiftY+fitsY*4+1]+
                               img[k,shiftX+fitsX*4 +2,shiftY+fitsY*4+2]+
                               img[k,shiftX+fitsX*4 +2,shiftY+fitsY*4+3]+
                               img[k,shiftX+fitsX*4 +3,shiftY+fitsY*4  ]+
                               img[k,shiftX+fitsX*4 +3,shiftY+fitsY*4+1]+
                               img[k,shiftX+fitsX*4 +3,shiftY+fitsY*4+2]+
                               img[k,shiftX+fitsX*4 +3,shiftY+fitsY*4+3])/16;
         img2[0,fitsX,fitsY]:=val/nrcolors;
    end;
  width2:=w;
  height2:=h;
  naxis3:=1;
end;


procedure bin_and_find_stars(img :image_array;binning:integer;cropping,hfd_min:double;get_hist{update hist}:boolean; out starlist3:star_list; out short_warning : string);{bin, measure background, find stars}
var
  old_width,old_height,old_naxis3,nrstars,i : integer;
  img_binned : image_array;

begin
  short_warning:='';{clear string}

  if ((binning>1) or (cropping<1)) then
  begin
    old_width:=width2;
    old_height:=height2;
    old_naxis3:=naxis3;
    if binning>1 then memo2_message('Creating grayscale x '+inttostr(binning)+' binning image for solving/star alignment.');
    if cropping<>1 then memo2_message('Cropping image x '+floattostrF2(cropping,0,2));

    if binning=2 then binX2_crop(cropping,img,img_binned) {combine values of 4 pixels, default option if 3 and 4 are not specified}
    else
    if binning=3 then binX3_crop(cropping,img,img_binned) {combine values of 9 pixels}
    else
    if binning=4 then binX4_crop(cropping,img,img_binned) {combine values of 16 pixels}
    else
    if binning=1 then binX1_crop(cropping,img,img_binned); {crop image, no binning}

    {test routine, to show bin result}
    //    img_loaded:=img_binned;
    //    naxis3:=1;
    //    plot_fits(mainwindow.image1,true);{plot real}
    //    exit;

    get_background(0,img_binned,true {load hist},true {calculate also standard deviation background},{var}cblack,star_level );{get back ground}
    find_stars(img_binned,hfd_min,starlist3); {find stars of the image and put them in a list}
    img_binned:=nil;
    nrstars:=Length(starlist3[0]);

    if height2<960 then
    begin
      short_warning:='Warning, remaining image dimensions too low! ';  {for FITS header and solution. Dimensions should be equal or better the about 1280x960}
      memo2_message('█ █ █ █ █ █ Warning, remaining image dimensions too low! Try to REDUCE OR REMOVE DOWNSAMPLING. Set this option in stack menu, tab alignment.');
    end;

    width2:=old_width; {restore to original size}
    height2:=old_height;
    naxis3:=old_naxis3;

    for i:=0 to nrstars-1 do {correct star positions for cropping. Simplest method}
    begin
      starlist3[0,i]:=starlist3[0,i]*binning+(width2*(1-cropping)/2);{correct star positions for binning/ cropping}
      starlist3[1,i]:=starlist3[1,i]*binning+(height2*(1-cropping)/2);
    end;
  end
  else
  begin
    if height2>2500 then
    begin
      short_warning:='Warning, increase downsampling!! '; {for FITS header and solution}
      memo2_message('Info: DOWNSAMPLING IS RECOMMENDED FOR LARGE IMAGES. Set this option in stack menu, tab alignment.');
    end
    else
    if height2<960 then
    begin
      short_warning:='Warning, small image dimensions! ';  {for FITS header and solution. Dimensions should be equal or better the about 1280x960}
      memo2_message('█ █ █ █ █ █ Warning, small image dimensions!');
    end;

    get_background(0,img,get_hist {load hist},true {calculate also standard deviation background}, {var} cblack,star_level);{get back ground}
    find_stars(img,hfd_min,starlist3); {find stars of the image and put them in a list}
  end;
end;

end.

