unit astap_stack;

{$mode ObjFPC}{$H+}

interface

uses
  Classes, SysUtils, math, fpjson,jsonparser,unit_command_line_general,unit_command_line_solving;



type
   theaderbackup  = record
     naxis   : integer;
     naxis1  : integer;
     naxis2  : integer;
     naxis3  : integer;
     crpix1 : double;{could be modified by crop}
     crpix2 : double;
     crval1 : double;
     crval2 : double;
     crota1 : double;{for 90 degrees rotate}
     crota2 : double;
     cdelt1 : double;
     cdelt2 : double;
     cd1_1  : double;
     cd1_2  : double;
     cd2_1  : double;
     cd2_2  : double;
     date_obs: string; {for accurate asteroid plotting after manual stack}
     header : string;
   end;
var
   header_backup : array of theaderbackup;{dynamic so memory can be freed}
  referenceX,referenceY    : double;{reference position used stacking}
  ref_X, ref_Y             : double;{reference position from FITS header, used for manual stacking of colour lights, second stage}

   {Variables for solving}
  SIN_dec0,COS_dec0,x_new_float,y_new_float,ra_ref,dec_ref,SIN_dec_ref,COS_dec_ref,crpix1_ref, crpix2_ref, CD1_1_ref, CD1_2_ref,CD2_1_ref,CD2_2_ref,exposure_ref,
  ap_0_1_ref,ap_0_2_ref,ap_0_3_ref,ap_1_0_ref,ap_1_1_ref, ap_1_2_ref,ap_2_0_ref,ap_2_1_ref,ap_3_0_ref, bp_0_1_ref,bp_0_2_ref,bp_0_3_ref,bp_1_0_ref,bp_1_1_ref,bp_1_2_ref,bp_2_0_ref,bp_2_1_ref,bp_3_0_ref   : double;

  a_order,ap_order: integer;{Simple Imaging Polynomial use by astrometry.net, if 2 then available}
  a_0_0,   a_0_1, a_0_2,  a_0_3,  a_1_0,  a_1_1,  a_1_2,  a_2_0,  a_2_1,  a_3_0 : double; {SIP, Simple Imaging Polynomial use by astrometry.net, Spitzer}
  b_0_0,   b_0_1, b_0_2,  b_0_3,  b_1_0,  b_1_1,  b_1_2,  b_2_0,  b_2_1,  b_3_0 : double; {SIP, Simple Imaging Polynomial use by astrometry.net, Spitzer}
  ap_0_0, ap_0_1,ap_0_2, ap_0_3, ap_1_0, ap_1_1, ap_1_2, ap_2_0, ap_2_1, ap_3_0 : double;{SIP, Simple Imaging Polynomial use by astrometry.net}
  bp_0_0, bp_0_1,bp_0_2, bp_0_3, bp_1_0, bp_1_1, bp_1_2, bp_2_0, bp_2_1, bp_3_0 : double;{SIP, Simple Imaging Polynomial use by astrometry.net}

type
{ TStack }
TStack = class(TObject)
 protected
   filename:string;
 public
   constructor Create(fn:string);
   destructor Destroy; override;
 private
   {Parse JSON and execute command whitin}
   function Execute:boolean;
   {Average darks files into a masterdark image}
   function CreateMasterDark(data: TJSONObject): boolean;
   {Average flat and bias images into a masterflat image}
   function CreateMasterFlat(data: TJSONObject): boolean;
   {Stack light images using optional flat and dark calibration images}
   function Stack(data: TJSONObject): boolean;
   {Actual sigma-clip stacking function}
   function stack_sigmaclip(oversize:integer; var files_to_process : array of String; out counter : integer;
                                    master_dark:image_array;
                                    master_flat:image_array;
                                    make_osc_color1:boolean;
                                    quad_tolerance:float;
                                    use_astrometry:boolean):boolean;
   {this routine works with mono files but makes coloured files mono, so less suitable for commercial cameras producing coloured raw lights}
   procedure average(mess:string; file_list : array of string; file_count:integer; var img2: image_array);
   {apply dark, flat if required, renew if different exposure or ccd temp}
   procedure apply_dark_and_flat(master_dark:image_array;master_flat:image_array;img:image_array);
   {helper functions}
   procedure backup_solution;{backup solution}
   function  restore_solution: boolean;{restore solution}
   procedure initialise_var1;{set variables correct}
   procedure initialise_var2;{set variables correct}
   procedure astrometric_to_vector;{convert astrometric solution to vector solution}
 end;
{main entry point}

function stack(filename:string):integer;


implementation

{ TStack }


{this routine works with mono files but makes coloured files mono, so less suitable for commercial cameras producing coloured raw lights}
procedure TStack.average(mess:string; file_list : array of string; file_count:integer; var img2: image_array);{combine to average or mean, make also mono from three colors if color}
var
   c,fitsX, fitsY : integer;
   img_tmp1 :image_array;
   img_loaded:boolean;
begin
  {average}
  for c:=0 to file_count-1 do
  begin
    writeln('Adding '+mess+' image '+inttostr(c+1)+ ' to '+mess+' average. '+file_list[c]);

    {load image}
    img_loaded:=load_fits(file_list[c],img_tmp1);

    if c=0 then {init}
    begin
      setlength(img2,1,width2,height2);{set length of image array mono}
      for fitsY:=0 to height2-1 do
        for fitsX:=0 to width2-1 do
         img2[0,fitsX,fitsY]:=0; {clear img}
    end;

    if naxis3=3 then  {for the rare case the darks are coloured. Should normally be not the case since it expects raw mono FITS files without bayer matrix applied !!}
    begin {color average}
      for fitsY:=0 to height2-1 do
         for fitsX:=0 to width2-1 do
           img2[0,fitsX,fitsY]:=img2[0,fitsX,fitsY]+(img_tmp1[0,fitsX,fitsY]+img_tmp1[1,fitsX,fitsY]+img_tmp1[2,fitsX,fitsY])/3;{fill with image}
    end
    else
    begin {mono average}
      for fitsY:=0 to height2-1 do
         for fitsX:=0 to width2-1 do
           img2[0,fitsX,fitsY]:=img2[0,fitsX,fitsY]+img_tmp1[0,fitsX,fitsY];{fill with image}

    end;
  end;{open files}

 if file_count>1 then {not required for single/master files}
  For fitsY:=0 to height2-1 do
     for fitsX:=0 to width2-1 do
       img2[0,fitsX,fitsY]:=img2[0,fitsX,fitsY]/file_count;{scale to one image}

  img_tmp1:=nil;{free mem}
end;


constructor TStack.Create(fn:string);
begin
filename:=fn;
end;

function TStack.Execute:boolean;
var
  f: TextFile;
  s: string;
  jsontext:string;
  jData : TJSONData;
  jObject : TJSONObject;
begin
assign(f,filename);
// Embed the file handling in a try/except block to handle errors gracefully
 try
   // Open the file for reading
   reset(f);

   // Keep reading lines until the end of the file is reached
   jsontext:='';
   while not eof(f) do
   begin
     readln(f, s);
     jsontext+=s;
   end;
   // Done so close the file
   CloseFile(f);
 except
   on E: EInOutError do
    writeln('File handling error occurred. Details: ', E.Message);
 end;
 jData := GetJSON(jsontext);
 jObject := TJSONObject(jData);
 // Dispatch different commands
 s := jObject.Get('Command');
 case (s) of
    'CreateMasterDark': self.CreateMasterDark(jObject);
    'CreateMasterFlat': self.CreateMasterFlat(jObject);
    'Stack': self.Stack(jObject);
 end;
end;

destructor TStack.Destroy;
begin
  inherited Destroy;
end;

function TStack.CreateMasterDark(data: TJSONObject):boolean;
var
  i:integer;
  fileCount:integer;
  file_list : array of string;
  Items: TJSonData;
  file_loaded: boolean;
  imgName:string;
  master_dark:image_array;
  master_dark_path:string;
begin
  Items := data.GetPath('Darks');
  fileCount:=Items.Count;
  file_list:=nil;
  setlength(file_list,fileCount);
  for i:=0 to fileCount-1 do
      file_list[i]:=data.FindPath(Format('Darks[%d]',[i])).AsString;
  master_dark:=nil;
  average('dark',file_list,fileCount,master_dark);

 // output file name
 if data.Get('Output')=Null then
   master_dark_path:=extractfilepath(file_list[0])+'master_dark_'+inttostr(fileCount)+'.fit'
 else master_dark_path:=data.Get('Output');


  if save_fits(master_dark,master_dark_path,-32,false) then
       begin
       writeln(Format('Success: Master dark image saved as %s',[master_dark_path]));
       Result:=True;
       end
  else begin
       writeln(Format('Error: Couldn''t save image at %s',[master_dark_path]));
       Result:=False;
       end;
  master_dark:=nil; {free image}
end;

function TStack.CreateMasterFlat(data: TJSONObject): boolean;
var
  fitsX,fitsY,i:integer;
  flatCount,biasCount:integer;
  file_list : array of string;
  Items: TJSonData;
  file_loaded: boolean;
  master_bias:image_array;
  master_flat:image_array;
  master_flat_path:string;
  biasWidth,biasHeight:integer;
  flatWidth,flatHeight:integer;
begin
  {Average DarkFlats (bias) if present}
  Items := data.GetPath('DarkFlats');
  biasCount:=Items.Count;
  if biasCount>0 then
      begin
      setlength(file_list,biasCount);
      for i:=0 to biasCount-1 do
          file_list[i]:=data.FindPath(Format('DarkFlats[%d]',[i])).AsString;
      average('darkflat',file_list,biasCount,master_bias);
      biasWidth:=length(master_bias[0]);{width}
      biasHeight:=length(master_bias[0,0]);{length}
      end
  else writeln('Warning: no flat-dark/bias found.');
  {Average Flats}
  Items := data.GetPath('Flats');
  flatCount:=Items.Count;
  if flatCount=0 then
      begin
      writeln('Warning: no flat files found.');
      abort;
      end;
  setlength(file_list,flatCount);
  for i:=0 to flatCount-1 do
      file_list[i]:=data.FindPath(Format('Flats[%d]',[i])).AsString;
  average('flat',file_list,flatCount,master_flat);
  flatWidth:=length(master_flat[0]);{width}
  flatHeight:=length(master_flat[0,0]);{length}
  { Apply bias if specified}
  if biasCount<>0 then
    begin
      if (flatWidth=biasWidth) and (flatHeight=biasHeight) then
        begin
        writeln('Using darkflat (bias) files on flat calibration.');
        for fitsY:=0 to height2-1 do
            for fitsX:=0 to width2-1 do
                begin
                {flats and bias already many mono in procedure average}
                master_flat[0,fitsX,fitsY]:=master_flat[0,fitsX,  fitsY  ] - master_bias[0,fitsX,  fitsY  ];
                end;
        end;
    end;

  {output file name}
  if data.Get('Output')=Null then
      master_flat_path:=extractfilepath(file_list[0])+'master_flat_'+inttostr(flatCount)+'.fit'
    else master_flat_path:=data.Get('Output');

  if save_fits(master_flat,master_flat_path,-32,false) then
       begin
       writeln(Format('Success: Master flat image saved as %s',[master_flat_path]));
       Result:=True;
       end
  else begin
       writeln(Format('Error: Couldn''t save image at %s',[master_flat_path]));
       Result:=False;
       end;
  master_flat:=nil; {free image}
  master_bias:=nil; {free image}
end;


{apply dark and flat if required, renew if different exposure or ccd temp}
procedure TStack.apply_dark_and_flat(master_dark:image_array;master_flat:image_array;img:image_array);
var
  fitsX,fitsY,k,light_naxis3, light_width,light_height,light_set_temperature : integer;
  calstat_local,
  light_date_obs              : string;
  datamax_light ,light_exposure,value,flat_factor,flat_norm_value,dark_norm_value,flat11,flat12,flat21,flat22 : double;
  img_loaded: boolean;
  darkWidth,darkHeight,flatWidth,flatHeight,imgWidth,imgHeight:integer;
begin
 {
  datamax_light:=datamax_org;

  light_naxis3:=naxis3; {preserve so it is not overriden by load dark_flat which will reset variable in load_fits}
  light_exposure:=exposure;{preserve so it is not overriden by apply dark_flat}
  light_width:=width2;
  light_height:=height2;
  light_date_obs:=date_obs;{preserve light date_obs}
  }



  imgWidth:=length(img[0]);{width}
  imgHeight:=length(img[0,0]);{length}

  height2:=imgHeight;
  width2:=imgWidth;


   if master_dark<>nil then begin
   darkWidth:=length(master_dark[0]);{width}
   darkHeight:=length(master_dark[0,0]);{length}
   if (imgWidth<>darkWidth) or (imgHeight<>darkHeight) then
      begin
      memo2_message('Error: image size mismatch!');
      memo2_message(Format('Master Dark: %d,%d pixels',[darkWidth,darkHeight]));
      memo2_message(Format('Image: %d,%d pixels',[imgWidth,imgHeight]));
      Halt;
      end;

   dark_norm_value:=0;
   for fitsY:=-2 to 3 do {do even times, 6x6}
       for fitsX:=-2 to 3 do
           dark_norm_value:=dark_norm_value+master_dark[0,fitsX+(width2 div 2),fitsY +(height2 div 2)];
      dark_norm_value:=round(dark_norm_value/36);{scale factor to apply flat. The norm value will result in a factor one for the center.}

   for fitsY:=0 to height2-1 do  {apply the dark}
      for fitsX:=0 to width2-1  do
      begin
        value:=master_dark[0,fitsX,fitsY]; {Darks are always made mono when making master dark}
        for k:=0 to naxis3-1 do {do all colors}
                    img[k,fitsX,fitsY]:=img[k,fitsX,fitsY] - value;
      end;
    calstat_local:=calstat_local+'D'; {dark applied}
    datamax_light:=datamax_light-dark_norm_value;
   end;

  if master_flat<>nil then
  begin
      flatWidth:=length(master_flat[0]);{width}
      flatHeight:=length(master_flat[0,0]);{length}
      if (imgWidth<>flatWidth) or (imgHeight<>flatHeight) then
         begin
         memo2_message('Error: image size mismatch!');
         memo2_message(Format('Master flat: %d,%d pixels',[flatWidth,flatHeight]));
         memo2_message(Format('Image: %d,%d pixels',[imgWidth,imgHeight]));
         Halt;
         end;

      flat_norm_value:=0;
      flat11:=0;
      flat12:=0;
      flat21:=0;
      flat22:=0;

      for fitsY:=-2 to 3 do {do even times, 6x6}
         for fitsX:=-2 to 3 do
         begin
           value:=master_flat[0,fitsX+(width2 div 2),fitsY +(height2 div 2)];
           flat_norm_value:=flat_norm_value+value;
           if ((odd(fitsX)) and (odd(fitsY)) ) then
                                 flat11:=flat11+value;
           if ((odd(fitsX)=false) and (odd(fitsY)) ) then
           flat12:=flat12+value;
           if ((odd(fitsX)) and (odd(fitsY)=false) ) then
           flat21:=flat21+value;
           if ((odd(fitsX)=false) and (odd(fitsY)=false) ) then
           flat22:=flat22+value;
         end;
      flat_norm_value:=round(flat_norm_value/36);{scale factor to apply flat. The norm value will result in a factor one for the center.}

      if max(max(flat11,flat12),max(flat21,flat22))/min(min(flat11,flat12),min(flat21,flat22))>2.0 then memo2_message('█ █ █ █ █ █ Warning flat pixels differ too much. Use white light for OSC flats or consider using option "Normalise OSC flat" █ █ █ █ █ █ ');

      for fitsY:=1 to height2 do  {apply the flat}
        for fitsX:=1 to width2 do
        begin
          flat_factor:=flat_norm_value/(master_flat[0,fitsX-1,fitsY-1]+0.001); {bias is already combined in flat in combine_flat}
          if abs(flat_factor)>3 then flat_factor:=1;{un-used sensor area? Prevent huge gain of areas only containing noise and no flat-light value resulting in very strong disturbing noise or high value if dark is missing. Typical problem for converted RAW's by Libraw}
          for k:=0 to naxis3-1 do {do all colors}
            img[k,fitsX-1,fitsY-1]:=img[k,fitsX-1,fitsY-1]*flat_factor;
        end;
    end;{flat correction}

end;


function TStack.Stack(data: TJSONObject): boolean;
var
  master_dark,master_flat:image_array;
  img:image_array;
  filename2,alignMethod: string;
  useAstrometry:boolean;
  img_loaded_ok:boolean;
  lightsCount,i:integer;
  file_list : array of string;
  Items: TJSonData;
begin
 // Get master dark
 if data.Get('MasterDark')<>Null then
    begin
    filename2:=data.Get('MasterDark');
     {load image}
    img_loaded_ok:=load_fits(filename2,master_dark);
    if img_loaded_ok=False then
       begin
       master_dark:=nil;
       memo2_message(Format('Warning: Couldn''t load %s!',[filename2]));
       end;
    end
 else  master_dark:=nil;
 // Get master flat
 if data.Get('MasterFlat')<>Null then
    begin
    filename2:=data.Get('MasterFlat');
     {load image}
    img_loaded_ok:=load_fits(filename2,master_flat);
    if img_loaded_ok=False then
       begin
       master_flat:=nil;
       memo2_message(Format('Warning: Couldn''t load %s!',[filename2]));
       end;
    end
 else  master_flat:=nil;
 // Get Lights
 Items := data.GetPath('Lights');
 lightsCount:=Items.Count;
 if lightsCount>0 then
     begin
     setlength(file_list,lightsCount);
     for i:=0 to lightsCount-1 do
         file_list[i]:=data.FindPath(Format('Lights[%d]',[i])).AsString;
     alignMethod := data.Get('Align');
     if alignMethod='stars' then useAstrometry:=False else useAstrometry:=True;
     if stack_sigmaclip(1, file_list,lightsCount,master_dark,master_flat,false,0.01,useAstrometry)=false then
        begin
          writeln('Error stacking images.');
          Result:=False;
          exit;
        end;
     {output file name}
     if data.Get('Output')=Null then
         filename2:=extractfilepath(file_list[0])+'stacked_'+inttostr(lightsCount)+'.fit'
       else filename2:=data.Get('Output');

     if save_fits(img_loaded,filename2,-32,false) then
          begin
          writeln(Format('Success: Stacked image saved as %s',[filename2]));
          Result:=True;
          end
     else begin
          writeln(Format('Error: Couldn''t save image at %s',[filename2]));
          Result:=False;
          end;
     end;

end;


procedure TStack.backup_solution;{backup solution}
begin
  if header_backup=nil then setlength(header_backup,1);{create memory}
  header_backup[0].crpix1:=crpix1;
  header_backup[0].crpix2:=crpix2;
  header_backup[0].crval1:=ra0;
  header_backup[0].crval2:=dec0;
  header_backup[0].crota1:=crota1;
  header_backup[0].crota2:=crota2;
  header_backup[0].cdelt1:=cdelt1;
  header_backup[0].cdelt2:=cdelt2;
  header_backup[0].cd1_1:=cd1_1;
  header_backup[0].cd1_2:=cd1_2;
  header_backup[0].cd2_1:=cd2_1;
  header_backup[0].cd2_2:=cd2_2;
  header_backup[0].date_obs:=date_obs;
//  if full then header_backup[0].header:=mainwindow.Memo1.Text;{backup fits header}
end;

function TStack.restore_solution: boolean;{restore solution and header memo}
begin
  if header_backup=nil then begin result:=false; exit;end;{no backup}

  crpix1:=header_backup[0].crpix1;
  crpix2:=header_backup[0].crpix2;

  ra0:=header_backup[0].crval1;
  dec0:=header_backup[0].crval2;

  crota1:=header_backup[0].crota1;
  crota2:=header_backup[0].crota2;
  cdelt1:=header_backup[0].cdelt1;
  cdelt2:=header_backup[0].cdelt2;
  cd1_1:=header_backup[0].cd1_1;
  cd1_2:=header_backup[0].cd1_2;
  cd2_1:=header_backup[0].cd2_1;
  cd2_2:=header_backup[0].cd2_2;
  date_obs:=header_backup[0].date_obs;
//  if full then mainwindow.Memo1.Text:=header_backup[0].header;{restore fits header}

  header_backup:=nil; {release memory}
  result:=true;
end;

procedure TStack.initialise_var1;{set variables correct}
begin
  ra_ref:=ra0;
  dec_ref:=dec0;
  sincos(dec_ref,SIN_dec_ref,COS_dec_ref);{do this in advance since it is for each pixel the same}
  crpix1_ref:=crpix1;
  crpix2_ref:=crpix2;
  CD1_1_ref:=CD1_1;
  CD1_2_ref:=CD1_2;
  CD2_1_ref:=CD2_1;
  CD2_2_ref:=CD2_2;

  exposure_ref:=exposure;
end;

procedure TStack.initialise_var2;{set variables correct}
begin
  ap_0_1_ref:=ap_0_1;{store polynomial first fits }
  ap_0_2_ref:=ap_0_2;
  ap_0_3_ref:=ap_0_3;
  ap_1_0_ref:=ap_1_0;
  ap_1_1_ref:=ap_1_1;
  ap_1_2_ref:=ap_1_2;
  ap_2_0_ref:=ap_2_0;
  ap_2_1_ref:=ap_2_1;
  ap_3_0_ref:=ap_3_0;

  bp_0_1_ref:=bp_0_1;
  bp_0_2_ref:=bp_0_2;
  bp_0_3_ref:=bp_0_3;
  bp_1_0_ref:=bp_1_0;
  bp_1_1_ref:=bp_1_1;
  bp_2_1_ref:=bp_2_1;
  bp_2_0_ref:=bp_2_0;
  bp_2_1_ref:=bp_2_1;
  bp_3_0_ref:=bp_3_0;
end;

procedure  calc_newx_newy(vector_based : boolean; fitsXfloat,fitsYfloat: double);  {apply either vector or astrometric correction}
var
  u,u0,v,v0,dRa,dDec,delta,ra_new,dec_new,delta_ra,det,gamma,SIN_dec_new,COS_dec_new,SIN_delta_ra,COS_delta_ra,h: double;
Begin

  if vector_based then {vector based correction}
  begin
     x_new_float:=solution_vectorX[0]*(fitsxfloat-1)+solution_vectorX[1]*(fitsYfloat-1)+solution_vectorX[2]; {correction x:=aX+bY+c  x_new_float in image array range 0..width2-1}
     y_new_float:=solution_vectorY[0]*(fitsxfloat-1)+solution_vectorY[1]*(fitsYfloat-1)+solution_vectorY[2]; {correction y:=aX+bY+c}
  end
  else
  begin {astrometric based correction}
    {6. Conversion (x,y) -> (RA,DEC)  for image to be added}
    u0:=fitsXfloat-crpix1;
    v0:=fitsYfloat-crpix2;

    if a_order>=2 then {apply SIP correction up second order}
    begin
      u:=u0 + a_2_0*u0*u0 + a_0_2*v0*v0 + a_1_1*u0*v0; {SIP correction}
      v:=v0 + b_2_0*u0*u0 + b_0_2*v0*v0 + b_1_1*u0*v0; {SIP correction}
    end
    else
    begin
      u:=u0;
      v:=v0;
    end;

    dRa :=(cd1_1 * u +cd1_2 * v)*pi/180;
    dDec:=(cd2_1 * u +cd2_2 * v)*pi/180;
    delta:=COS_dec0 - dDec*SIN_dec0;
    gamma:=sqrt(dRa*dRa+delta*delta);
    RA_new:=ra0+arctan(Dra/delta);
    DEC_new:=arctan((SIN_dec0+dDec*COS_dec0)/gamma);


   {5. Conversion (RA,DEC) -> (x,y) of reference image}
    sincos(dec_new,SIN_dec_new,COS_dec_new);{sincos is faster then separate sin and cos functions}
    delta_ra:=RA_new-ra_ref;
    sincos(delta_ra,SIN_delta_ra,COS_delta_ra);

    H := SIN_dec_new*sin_dec_ref + COS_dec_new*COS_dec_ref*COS_delta_ra;
    dRA := (COS_dec_new*SIN_delta_ra / H)*180/pi;
    dDEC:= ((SIN_dec_new*COS_dec_ref - COS_dec_new*SIN_dec_ref*COS_delta_ra ) / H)*180/pi;

    det:=CD2_2_ref*CD1_1_ref - CD1_2_ref*CD2_1_ref;

    u0:= - (CD1_2_ref*dDEC - CD2_2_ref*dRA) / det;
    v0:= + (CD1_1_ref*dDEC - CD2_1_ref*dRA) / det;

    if ap_order>=2 then {apply SIP correction up to second order}
    begin
      x_new_float:=(CRPIX1_ref + u0+ap_0_1*v0+ ap_0_2*v0*v0+ ap_0_3*v0*v0*v0 +ap_1_0*u0 + ap_1_1*u0*v0+  ap_1_2*u0*v0*v0+ ap_2_0*u0*u0 + ap_2_1*u0*u0*v0+  ap_3_0*u0*u0*u0)-1;{3th order SIP correction, fits count from 1, image from zero therefore subtract 1}
      y_new_float:=(CRPIX2_ref + v0+bp_0_1*v0+ bp_0_2*v0*v0+ bp_0_3*v0*v0*v0 +bp_1_0*u0 + bp_1_1*u0*v0+  bp_1_2*u0*v0*v0+ bp_2_0*u0*u0 + bp_2_1*u0*u0*v0+  bp_3_0*u0*u0*u0)-1;{3th order SIP correction}
    end
    else
    begin
      x_new_float:=(CRPIX1_ref + u0)-1; {in image array range 0..width-1}
      y_new_float:=(CRPIX2_ref + v0)-1;
    end;
  end;{astrometric}
end;{calc_newx_newy}

procedure TStack.astrometric_to_vector;{convert astrometric solution to vector solution}
var
  flipped,flipped_reference  : boolean;

begin
  a_order:=0; {SIP correction should be zero by definition}

  calc_newx_newy(false,1,1) ;
  solution_vectorX[2]:=x_new_float;
  solution_vectorY[2]:=y_new_float;

  calc_newx_newy(false,2, 1); {move one pixel in X}

  solution_vectorX[0]:=+(x_new_float- solution_vectorX[2]);
  solution_vectorX[1]:=-(y_new_float- solution_vectorY[2]);

  calc_newx_newy(false,1, 2);{move one pixel in Y}
  solution_vectorY[0]:=-(x_new_float- solution_vectorX[2]);
  solution_vectorY[1]:=+(y_new_float- solution_vectorY[2]);

  if (cd2_2)<>0 then
    begin
    flipped:=(cd1_1/cd2_2)>0; {flipped image. Either flipped vertical or horizontal but not both. Flipped both horizontal and vertical is equal to 180 degrees rotation and is not seen as flipped}
    flipped_reference:=(cd1_1_ref/cd2_2_ref)>0; {flipped reference image}
    if flipped<>flipped_reference then {this can happen is user try to add images from a diffent camera/setup}
         begin
         solution_vectorX[1]:=-solution_vectorX[1];
         solution_vectorY[0]:=-solution_vectorY[0];
         end;
    end;
end;

{stack using sigma clip average}
function TStack.stack_sigmaclip(oversize:integer; var files_to_process : array of String; out counter : integer;
                                 master_dark:image_array;
                                 master_flat:image_array;
                                 make_osc_color1:boolean;
                                 quad_tolerance:float;
                                 use_astrometry:boolean):boolean;
type
   tsolution  = record
     solution_vectorX : solution_vector {array[0..2] of double};
     solution_vectorY : solution_vector;
     cblack : double;
   end;
var
    solutions      : array of tsolution;
    fitsX,fitsY,c,width_max, height_max, old_width, old_height,x_new,y_new,col ,binning,oversizeV         : integer;
    background_correction, variance_factor, value,weightF,hfd_min                                         : double;
    init, solution:boolean;
    use_star_alignment, use_astrometry_internal:boolean;
    vector_based :boolean;
    warning  : string;
    img_loaded_ok : boolean;
    pedestal_s : double;{target background value}
begin
    {move often uses setting to booleans. Great speed improved if use in a loop and read many times}
    variance_factor:=2.0;//sqr(strtofloat2(stackmenu1.sd_factor1.text));

    {to ignore hot pixels which are too small}

    hfd_min:=0.1; //max(0.8 {two pixels},strtofloat2(stackmenu1.min_star_size_stacking1.caption){hfd});

    {Select alignment method}
    use_star_alignment:=False;
    use_astrometry_internal:=False;
    if   use_astrometry=True then
         use_astrometry_internal:=True
    else use_star_alignment:=True;

    counter:=0;

    init:=false;
    background_correction:=0;{required for astrometric alignment}
    {light average}
    begin
      counter:=0;
      setlength(solutions,length(files_to_process));
      init:=false;
      for c:=0 to length(files_to_process)-1 do
      if length(files_to_process[c])>0 then
      begin
      try { Do some lengthy operation }
        filename2:=files_to_process[c];
        writeln(Format('Sigma clip stacking: Processing file %s',[filename2]));

        {load image}
        img_loaded_ok:=load_fits(filename2,img_loaded);
        if img_loaded_ok=false then begin memo2_message('Abort!! Can'+#39+'t load '+filename2); exit(false);end;

        {Solve image if required by astrometry align}
        if use_astrometry_internal=true then
          if cd1_1=0 then
             begin
             memo2_message('Image '+filename2+' not solved. Solving now...');
             solve_image(img_loaded,true);
             SaveFITSwithupdatedheader1;
             end;

        width2:=length(img_loaded);{width}
        height2:=length(img_loaded);{length}
        naxis3:=length(img_loaded); {colors}

        if init=true then
        begin
           if ((old_width<>width2) or (old_height<>height2)) then memo2_message('Warning: different size image!');
           if naxis3>length(img_average) {naxis3} then begin memo2_message('Abort!! Can'+#39+'t combine mono and colour files.'); exit(false);end;
        end;


        if init=false then {phase (1) average}
        begin
          binning:=report_binning(height2);{select binning based on the height of the light}
          backup_solution;{backup solution}
          initialise_var1;{set variables correct}
          initialise_var2;{set variables correct}

//          if ((bayerpat='') and (make_osc_color1.checked)) then
//             if stackmenu1.bayer_pattern1.Text='auto' then memo2_message('Warning, Bayer colour pattern not in the header! Check colours and if wrong set Bayer pattern manually in tab "stack alignment". █ █ █ █ █ █')
//             else
//             if test_bayer_matrix(img_loaded)=false then  memo2_message('Warning, grayscale image converted to colour! Un-check option "convert OSC to colour". █ █ █ █ █ █');
        end;

        {apply dark, flat if required, renew if different exposure or ccd temp}
        apply_dark_and_flat(master_dark,master_flat,img_loaded);

        {these global variables are passed-on in procedure to protect against overwriting}

        memo2_message('Adding light file: '+filename2+' dark compensated to light average.');

        if make_osc_color1 then {do demosaic bayer}
        begin
          if naxis3>1 then memo2_message('Warning: light is already in colour ! Will skip demosaic. █ █ █ █ █ █')
          else
          //demosaic_bayer(img_loaded); {convert OSC image to colour} TODO
         {naxis3 is now 3}
        end;

        if ((init=false ) and (use_astrometry_internal=false)) then {first image and not astrometry_internal}
        begin
            bin_and_find_stars(img_loaded, binning,1  {cropping},hfd_min,true{update hist},starlist1,warning);{bin, measure background, find stars}

            find_quads(starlist1,0,quad_smallest,quad_star_distances1);{find quads for reference image}
            pedestal_s:=cblack;{correct for difference in background, use cblack from first image as reference. Some images have very high background values up to 32000 with 6000 noise, so fixed pedestal_s of 1000 is not possible}
            if pedestal_s<500 then pedestal_s:=500;{prevent image noise could go below zero}
            background_correction:=pedestal_s-cblack;
            datamax_org:=datamax_org+background_correction; if datamax_org>$FFFF then  datamax_org:=$FFFF; {note datamax_org is already corrected in apply dark}
        end;

        if init=false then {init}
        begin
          if oversize<0 then {shrink, adapt in ratio}
          begin
            oversize:=max(oversize,-round((width2-100)/2) );{minimum image width is 100}
            oversizeV:=round(oversize*height2/width2);{vertical}
            height_max:=height2+oversizeV*2;
          end
          else
          begin
            oversizeV:=oversize;
            height_max:=height2+oversize*2;
          end;
          width_max:=width2+oversize*2;

          setlength(img_average,naxis3,width_max,height_max);
          setlength(img_temp,naxis3,width_max,height_max);
            for fitsY:=0 to height_max-1 do
             for fitsX:=0 to width_max-1 do
               for col:=0 to naxis3-1 do
               begin
                 img_average[col,fitsX,fitsY]:=0; {clear img_average}
                 img_temp[col,fitsX,fitsY]:=0; {clear img_temp}
               end;
          old_width:=width2;
          old_height:=height2;
        end;{init, c=0}

        solution:=true;
        if use_astrometry_internal then sincos(dec0,SIN_dec0,COS_dec0) {do this in advance since it is for each pixel the same}
        else
        begin {align using star match}
           if init=true then {second image}
                begin{internal alignment}
                  bin_and_find_stars(img_loaded, binning,1  {cropping},hfd_min,true{update hist},starlist2,warning);{bin, measure background, find stars}

                  background_correction:=pedestal_s-cblack;{correct later for difference in background}
                  datamax_org:=datamax_org+background_correction; if datamax_org>$FFFF then  datamax_org:=$FFFF; {note datamax_org is already corrected in apply dark}

                  find_quads(starlist2,0,quad_smallest,quad_star_distances2);{find star quads for new image}
                  if find_offset_and_rotation(3,quad_tolerance) then {find difference between ref image and new image}
                  begin
                    memo2_message(inttostr(nr_references)+' of '+ inttostr(nr_references2)+' quads selected matching within '+FloatToStr(quad_tolerance)+' tolerance.'
                       +'  Solution x:='+floattostr6(solution_vectorX[0])+'*x+ '+floattostr6(solution_vectorX[1])+'*y+ '+floattostr6(solution_vectorX[2])
                       +',  y:='+floattostr6(solution_vectorY[0])+'*x+ '+floattostr6(solution_vectorY[1])+'*y+ '+floattostr6(solution_vectorY[2]) );

                    solutions[c].solution_vectorX:= solution_vectorX;{store solutions}
                    solutions[c].solution_vectorY:= solution_vectorY;
                    solutions[c].cblack:=cblack;


                  end
                    else
                    begin
                      memo2_message('Not enough quad matches <3 or inconsistent solution, skipping this image.');
                      solution:=false;
                    end;
                 end{internal alignment}
              else
              begin {first image}
                reset_solution_vectors(1);{no influence on the first image}
                solutions[c].solution_vectorX:= solution_vectorX; {store solutions for later}
                solutions[c].solution_vectorY:= solution_vectorY;
                solutions[c].cblack:=cblack;
               end;

        end;
        init:=true;{initialize for first image done}

        if solution then
        begin

          inc(counter);
          if exposure<>0 then weightF:=exposure/exposure_ref else weightF:=1;{influence of each image depending on the exposure_time}

          vector_based:=(use_star_alignment);
          if ((vector_based=false) and (a_order=0)) then {no SIP from astronomy.net}
          begin
            astrometric_to_vector;{convert astrometric solution to vector solution}
            vector_based:=true;
          end;

          for fitsY:=1 to height2 do {skip outside "bad" pixels if mosaic mode}
          for fitsX:=1 to width2  do
          begin
            calc_newx_newy(vector_based,fitsX,fitsY);{apply correction}
            x_new:=round(x_new_float+oversize);y_new:=round(y_new_float+oversizeV);

            if ((x_new>=0) and (x_new<=width_max-1) and (y_new>=0) and (y_new<=height_max-1)) then
            begin
              for col:=0 to naxis3-1 do
              begin
                img_average[col,x_new,y_new]:=img_average[col,x_new,y_new]+ img_loaded[col,fitsX-1,fitsY-1]+background_correction;{Note fits count from 1, image from zero}
                img_temp[col,x_new,y_new]:=img_temp[col,x_new,y_new]+weightF {norm 1};{count the number of image pixels added=samples}
              end;
            end;
          end;
        end;
        finally
        end;
      end;{try}
      if counter<>0 then
      For fitsY:=0 to height_max-1 do
        for fitsX:=0 to width_max-1 do
            for col:=0 to naxis3-1 do
            if img_temp[col,fitsX,fitsY]<>0 then
               img_average[col,fitsX,fitsY]:=img_average[col,fitsX,fitsY]/img_temp[col,fitsX,fitsY];{scale to one image by diving by the number of pixels added}

    end;  {light average}

    {standard deviation of light images}  {stack using sigma clip average}
    begin {standard deviation}
      counter:=0;

      init:=false;
      for c:=0 to length(files_to_process)-1 do
      if length(files_to_process[c])>0 then
      begin
        try { Do some lengthy operation }
          filename2:=files_to_process[c];

          {load image}
          img_loaded_ok:=load_fits(filename2,img_loaded);
          width2:=length(img_loaded);{width}
          height2:=length(img_loaded);{length}
          naxis3:=length(img_loaded); {colors}

          if init=false then
          begin
            {not required. Done in first step}
          end;

          {apply dark, flat if required, renew if different exposure or ccd temp}
          apply_dark_and_flat(master_dark,master_flat,img_loaded);

          memo2_message('Calculating pixels σ of light file '+inttostr(c+1)+'-'+filename2);

          {if make_osc_color1.checked then {do demosaic bayer}
          begin
            if naxis3>1 then memo2_message('█ █ █ █ █ █ Warning, light is already in colour ! Will skip demosaic. █ █ █ █ █ █')
            else
            demosaic_bayer(img_loaded); {convert OSC image to colour}
           {naxis3 is now 3}
          end;}

          if init=false then {init (2) for standard deviation step}
          begin
            setlength(img_variance,naxis3,width_max,height_max);{mono}
            for fitsY:=0 to height_max-1 do
            for fitsX:=0 to width_max-1 do
            begin
              for col:=0 to naxis3-1 do img_variance[col,fitsX,fitsY]:=0; {clear img_average}
            end;
            old_width:=width2;
            old_height:=height2;
          end;{c=0}

          inc(counter);

          if use_astrometry_internal then  sincos(dec0,SIN_dec0,COS_dec0) {do this in advance since it is for each pixel the same}
          else
          begin {align using star match, read saved solution vectors}
              solution_vectorX:=solutions[c].solution_vectorX; {restore solution}
              solution_vectorY:=solutions[c].solution_vectorY;
              cblack:=solutions[c].cblack;
              background_correction:=pedestal_s-cblack;{correction for difference in background}
              datamax_org:=datamax_org+background_correction; if datamax_org>$FFFF then  datamax_org:=$FFFF; {note datamax_org is already corrected in apply dark}
          end;
          init:=true;{initialize for first image done}

          vector_based:=use_star_alignment;
          if ((vector_based=false) and (a_order=0)) then {no SIP from astronomy.net}
          begin
            astrometric_to_vector;{convert astrometric solution to vector solution}
            vector_based:=true;
          end;

          for fitsY:=1 to height2 do {skip outside "bad" pixels if mosaic mode}
          for fitsX:=1 to width2  do
          begin
            calc_newx_newy(vector_based,fitsX,fitsY);{apply correction}
            x_new:=round(x_new_float+oversize);y_new:=round(y_new_float+oversizeV);
            if ((x_new>=0) and (x_new<=width_max-1) and (y_new>=0) and (y_new<=height_max-1)) then
            begin
              for col:=0 to naxis3-1 do img_variance[col,x_new,y_new]:=img_variance[col,x_new,y_new] +  sqr( img_loaded[col,fitsX-1,fitsY-1]+ background_correction - img_average[col,x_new,y_new]); {Without flats, sd in sqr, work with sqr factors to avoid sqrt functions for speed}
            end;
          end;
        finally
        end;
      end;{try}
      if counter<>0 then
        For fitsY:=0 to height_max-1 do
          for fitsX:=0 to width_max-1 do
            for col:=0 to naxis3-1 do
              if img_temp[col,fitsX,fitsY]<>0 then {reuse the img_temp from light average}
                 img_variance[col,fitsX,fitsY]:=1+img_variance[col,fitsX,fitsY]/img_temp[col,fitsX,fitsY]; {the extra 1 is for saturated images giving a SD=0}{scale to one image by diving by the number of pixels tested}
    end; {standard deviation of light images}


    {throw out the outliers of light-dark images}  {stack using sigma clip average}
    begin
      counter:=0;
      init:=false;
      for c:=0 to length(files_to_process)-1 do
      if length(files_to_process[c])>0 then
      begin
        try { Do some lengthy operation }
          filename2:=files_to_process[c];

          {load image}
          img_loaded_ok:=load_fits(filename2,img_loaded);
          width2:=length(img_loaded);{width}
          height2:=length(img_loaded);{length}
          naxis3:=length(img_loaded); {colors}

          {apply dark, flat if required, renew if different exposure or ccd temp}
          apply_dark_and_flat(master_dark,master_flat,img_loaded);

          memo2_message('Combining '+inttostr(c+1)+'-'+filename2+'", ignoring outliers.');

          {if make_osc_color1.checked then {do demosaic bayer}
          begin
            if naxis3>1 then memo2_message('█ █ █ █ █ █ Warning, light is already in colour ! Will skip demosaic. █ █ █ █ █ █')
            else
            demosaic_bayer(img_loaded); {convert OSC image to colour}
           {naxis3 is now 3}
          end;}

          if init=false then {init, (3) step throw outliers out}
          begin
            setlength(img_temp,naxis3,width_max,height_max);
            setlength(img_final,naxis3,width_max,height_max);
            for fitsY:=0 to height_max-1 do
            for fitsX:=0 to width_max-1 do
            begin
              for col:=0 to naxis3-1 do
              begin
                img_temp[col,fitsX,fitsY]:=0; {clear img_temp}
                img_final[col,fitsX,fitsY]:=0; {clear img_temp}
              end;
            end;
            old_width:=width2;
            old_height:=height2;
          end;{init}

          inc(counter);

          if use_astrometry_internal then  sincos(dec0,SIN_dec0,COS_dec0) {do this in advance since it is for each pixel the same}
          else
          begin {align using star match, read saved solution vectors}
              solution_vectorX:=solutions[c].solution_vectorX; {restore solution}
              solution_vectorY:=solutions[c].solution_vectorY;
              cblack:=solutions[c].cblack;
              background_correction:=pedestal_s-cblack;{correct for difference in background}
              datamax_org:=datamax_org+background_correction; if datamax_org>$FFFF then  datamax_org:=$FFFF; {note datamax_org is already corrected in apply dark}
          end;
          init:=true;{initialize for first image done}

          vector_based:=use_star_alignment;
          if ((vector_based=false) and (a_order=0)) then {no SIP from astronomy.net}
          begin
            astrometric_to_vector;{convert astrometric solution to vector solution}
            vector_based:=true;
          end;

          for fitsY:=1 to height2 do
          for fitsX:=1 to width2  do
          begin
            calc_newx_newy(vector_based,fitsX,fitsY);{apply correction}
            x_new:=round(x_new_float+oversize);y_new:=round(y_new_float+oversizeV);
            if ((x_new>=0) and (x_new<=width_max-1) and (y_new>=0) and (y_new<=height_max-1)) then
            begin
              for col:=0 to naxis3-1 do {do all colors}
              begin
                value:=img_loaded[col,fitsX-1,fitsY-1]+ background_correction;
                if sqr (value - img_average[col,x_new,y_new])< variance_factor*{sd sqr}( img_variance[col,x_new,y_new])  then {not an outlier}
                begin
                  img_final[col,x_new,y_new]:=img_final[col,x_new,y_new]+ value;{dark and flat, flat dark already applied}
                  img_temp[col,x_new,y_new]:=img_temp[col,x_new,y_new]+weightF {norm 1};{count the number of image pixels added=samples}
                end;
              end;
            end;
          end;

          finally
        end;
      end;

     {scale to number of pixels}
      if counter<>0 then
      begin
        height2:=height_max;
        width2:=width_max;
        setlength(img_loaded,naxis3,width2,height2);{new size}

        for col:=0 to naxis3-1 do {do one or three colors} {compensate for number of pixel values added per position}
          For fitsY:=0 to height2-1 do
            for fitsX:=0 to width2-1 do
            if img_temp[col,fitsX,fitsY]<>0 then img_loaded[col,fitsX,fitsY]:=img_final[col,fitsX,fitsY]/img_temp[col,fitsX,fitsY] {scale to one image by diving by the number of pixels added}
            else
            begin { black spot filter. Note for this version img_temp is counting for each color since they could be different}
              if ((fitsX>0) and (fitsY>0)) then {black spot filter, fix black spots which show up if one image is rotated}
              begin
                if ((img_temp[col,fitsX-1,fitsY]<>0){and (img_temp[col,fitsX,fitsY-1]<>0)}{keep borders nice for last pixel right}) then img_loaded[col,fitsX,fitsY]:=img_loaded[col,fitsX-1,fitsY]{take nearest pixel x-1 as replacement}
                else
                if img_temp[col,fitsX,fitsY-1]<>0 then img_loaded[col,fitsX,fitsY]:=img_loaded[col,fitsX,fitsY-1]{take nearest pixel y-1 as replacement}
                else
                img_loaded[col,fitsX,fitsY]:=0;{clear img_loaded since it is resized}
              end {fill black spots}
              else
              img_loaded[col,fitsX,fitsY]:=0;{clear img_loaded since it is resized}
            end; {black spot filter}
      end;{counter<>0}

      restore_solution;{restore solution variable of reference image for annotation and mount pointer}

    end;{throw out the outliers of light-dark images}
  {image arrays will be nilled later. This is done for early exits}

  solutions:=nil;
  Result:=true;
end;   {stack using sigma clip average}



function stack(filename:string):integer;
var
  stacker:TStack;
begin
stacker:= TStack.Create(filename);
stacker.Execute;
Result:=0;
end;

end.


