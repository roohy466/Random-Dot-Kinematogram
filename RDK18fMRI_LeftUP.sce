
#screen_width=1280;
#screen_height=1024;
#screen_bit_depth = 32;
scenario = "RandomDotKinematogram";

no_logfile=false; 
response_matching = simple_matching;
active_buttons=2;                   ## HID keyboard device number 4 (left) and 7 (right) 
button_codes=1,2;

write_codes = true;  ##the output port is LPT3
pulse_width = 5;
response_logging = log_active ; 

default_font_size = 60;
default_background_color=0,0,0;
default_text_color=150,30,190;   



# ------ FMRI stuff --------
#scenario_type=fMRI_emulation;
scenario_type=fMRI;
scan_period=2200; 
pulses_per_scan=1;  
sequence_interrupt = true;  # CHECK THIS    
pulse_code=9;


begin;

# Fixation dot
ellipse_graphic{
   ellipse_width = 3;
   ellipse_height = 3;
   color = 255, 0, 0, 200;
}EllipseShape;

ellipse_graphic{
   ellipse_width = 3;
   ellipse_height = 3;
   color = 55, 255, 0, 200;
}EllipseShapeGreen;

ellipse_graphic{
   ellipse_width = 3;
   ellipse_height = 3;
   color = 0, 0, 0, 255;
}EllipseShapeBlack;

trial {
	stimulus_event { 
		picture{	
			ellipse_graphic EllipseShape;
			center_x = 0; 
			center_y = 0;
		} pic;
		code=11;
		port_code=11;
	} main_event;
} main_trial;

trial{ 
	trial_duration=3000;
      picture {
			text {caption=" "; font_color=255,255,255; font_size = 35;} timing; 
				x = 0; y = 70;
			text {caption=" "; font_color=255,255,255; font_size = 24;} timing1;
				x = 0; y = -150;
		} rest;
   time=0;
   code= "555";
   port_code=99;
}long_PAUSE;  #rest between two runs


trial { 
	#trial_type = first_response;         
	stimulus_event { 
		nothing{};
		time = 0;
      code= "666";
      port_code=77;
	} StimEvent;
} TrialName;#  write dummy trial for the event code in logfile 
trial {    
	stimulus_event { 
		nothing{};
		time = 0;
		code= "777";
	} DummyEvent;
} DummyTrial;#  write dummy trial for the event code in logfile 

###################          PCL          ###########################

begin_pcl;

# modality to measure
int fMRI_exp=1,   # for behavioral data is 0;
	ExpStartVolTime=0,
	padding_volumes=5, #padding time
	VolumeTR=2200; #msecs

# stim type
int up=1,   ######>>>>>>>>>>>>>>>>>>>>>>>>>>>.  defined the up down button  1 = left, 2 = right
	 down=2;


# Block,conditions, trial information
int num_trials_percond =8 ,#5; # number of trials for each blocks
	num_blocks=4;#16;  # number of blocks
#array<double>Coherence[1]={ 56.0}; #Coherence of each trial
array<double>Coherence[3]={32.0, 21.0, 15.0}; #Coherence of each trial
array<int>Orient[2] ={-1 ,1};#the orientation of each trial left=-1 , right=1 
#array<int>MSeq[48]={1,1,5,2,1,2,6,0,3,3,1,6,3,6,4,0,2,2,3,4,2,4,5,0,6,6,2,5,6,5,1,0,4,4,6,1,4,1,3,0,5,5,4,3,5,3,2,0};  # the M sequence for 49 trials
int SequenceNumber=168;
array<int>MSeqPerm[SequenceNumber]={1,1,12,10,12,5,7,10,9,2,10,6,12,0,2,2,11,7,11,10,1,7,5,4,7,12,11,0,4,4,9,1,9,7,2,1,10,8,1,11,9,0,
				8,8,5,2,5,1,4,2,7,3,2,9,5,0,3,3,10,4,10,2,8,4,1,6,4,5,10,0,6,6,7,8,7,4,3,8,2,12,8,10,7,0,12,12,1,3,1,8,6,
				3,4,11,3,7,1,0,11,11,2,6,2,3,12,6,8,9,6,1,2,0,9,9,4,12,4,6,11,12,3,5,12,2,4,0,5,5,8,11,8,12,9,11,6,10,11,
				4,8,0,10,10,3,9,3,11,5,9,12,7,9,8,3,0,7,7,6,5,6,9,10,5,11,1,5,3,6,0};# the M sequence for 169 trials
  # the M sequence for 49 trials
array<int>TMPMSeqPerm[SequenceNumber];
array<int>MSeq[SequenceNumber];

TMPMSeqPerm.fill(1,0,0,0);
int shiftsequence=random(1,SequenceNumber);
term.print_line(shiftsequence);

loop int i=1 until i > SequenceNumber   # number of orientations (here is two)
begin;
	if i<=shiftsequence then
		TMPMSeqPerm[(SequenceNumber-shiftsequence+i)]=MSeqPerm[i];
		#term.print_line("There " + string(SequenceNumber-shiftsequence+i)+ "   " +string(TMPMSeqPerm[((SequenceNumber-shiftsequence+i))]));
	else
		TMPMSeqPerm[(i-shiftsequence)]=MSeqPerm[i];
		#term.print_line("here " + string(i-shiftsequence)+ "   " +string(TMPMSeqPerm[(i-shiftsequence)]));
	end;	
	i=i+1;
end;
MSeqPerm.assign(TMPMSeqPerm);
loop int i=1 until i > SequenceNumber   # number of orientations (here is two)
begin;
	if MSeqPerm[i] > 6 then
		MSeq[i]=MSeqPerm[i]-6;
	else
		MSeq[i]=MSeqPerm[i];
	end;
	#term.print_line(MSeq[i]);
	i=i+1;
end;
int num_trials_perBlock =MSeq.count()/num_blocks;#5; # number of trials for each blocks
	
#array<int>PerfomanceDur[48]=
#array<int>MSeq[48]={6,4,0,2,2,3,4,2,4,5,0,6,6,2,5,6,5,1,0,4,4,6,1,4,1,3,0,5,5,4,3,5,3,2,0,1,1,5,2,1,2,6,0,3,3,1,6,3};  # the M sequence for 49 trials
#array<int>MSeq[48]={2,3,4,2,4,5,0,6,6,2,5,6,5,1,0,4,4,6,1,4,1,3,0,5,5,4,3,5,3,2,0,1,1,5,2,1,2,6,0,3,3,1,6,3,6,4,0,2};  
int num_cond=Coherence.count()*Orient.count(),
    num_trials = (num_cond+1)*num_trials_percond; 
#term.print_line(num_trials);

# Size of stimulus
double scrn_widthDegree =  6.7, # width of screen in degree
	# screen information and distance   
	screenComercialSize=27.0, # the comercial size in inches
	screen_Xresolution=1280.0, #in pixels
	screen_Yresolution=800.0, #in Pixels
	screen_Zdistance=200.0, # in centemiter
	XScr_size_in=screenComercialSize*cos( arctan(screen_Yresolution/screen_Xresolution)),
	# calculating the resolution of the screen
	SCRVD=(arctan(( 2.54* XScr_size_in / 2.0 ) / screen_Zdistance ) * 2.0)*180.0/pi_value,
	SCPixperVD=screen_Xresolution/SCRVD,
   scrn_width =  scrn_widthDegree*SCPixperVD, # width of screen in pixels
	
	Rmax= floor(sqrt(2.0*scrn_width*scrn_width)/2.1)- scrn_width/10.0;

# Random dots information
int num_dots= 600, #400;  # number of dots
	initial_random_frames=2000;  # Number of random location frames for all dots
double speed_degree = 8.0, # speed in degree
	SizeScale=screen_Zdistance/70.0, # 75 cm is for normal psychophysics
	FixEllipse = 18.0,  # R of Fixation Ellipse 
	dot_radius = 1.3*SizeScale,  # radius of each dot, which creates by number of mini lines = CircleDivision value
	line_width=1.8*SizeScale; # the width of mini line
int CircleDivision=2, # the number of mini line
	MaxLifeTime = 30,  #Maximum life time of a dot
	trial_duration=1500, # tral duration for memory reason
	accepted_trial_duration=1800; # Trial Duration


# Stmulus Timing 
int wait_interval_time_btwenn_frame= int(ceil((1000.0/60.0)*1.0)); # time interval between frames
int trial_duration_frames=int(floor(double(trial_duration/wait_interval_time_btwenn_frame)));
double move_amount = (speed_degree *double(wait_interval_time_btwenn_frame)*scrn_width) /(scrn_widthDegree*1000.0) ; # amount, in pixels, each dot moves on each vertial refresh (velocity)
term.print_line("Frames=  " + string( trial_duration_frames) +"  InterFrame=  " + string(wait_interval_time_btwenn_frame) + "  TotalDur=  " +string (wait_interval_time_btwenn_frame*trial_duration_frames));
array<int>InterTrialDealy[MSeq.count()];

# Random dot shape and initilization
int KDivision= int(floor(sqrt(double(num_dots/2)))); # the devision of R tothis number (K)value for Random dot position function
if KDivision>=8  then
	KDivision=8;  # Maximun devision can reach accordng to the size, the more slows down the calculation however it creates more homogeneous random dot picture
end; 
double CurrentAngle=0.0;
double CurrentCoherent=0.0;

array<double>Orient_Coherence[0][0]; # the coherence and orientation of each trials
array<int>lifetime_matrix[num_dots];  # this variable availabe in all subroutines
array<int>lifetime_matrix_fixed[num_dots]; 
array<double>TotalLogData[num_trials*num_blocks][6];  # LOG of whole experiment
array<double>rand_locations_frame[trial_duration_frames][num_dots][2]; # Location for all dots in the frame


array<double>TrailMotionEnergy_frame[trial_duration_frames][8]; # Motion energy in 8 quadrants 0 -> Rmax/2,  Rmax/2 -> Rmax
array<double>TrailMotionEnergy_trails[num_trials*num_blocks][trial_duration_frames][8];


#term.print_line(num_trials);
response_data LastResp;
stimulus_data LastStim,LastStimFirst; 

array<int>nums_permutation[num_dots]; #counting the number a number repeated
nums_permutation.fill(1,0,1,1);

# orientation of movining dot , Random Angels
int CheckFirstTrial=0; 
int RangDevision= num_trials_percond*num_blocks;

double fillArrayJitter= 30.0/double(RangDevision); # the slice covering for top or down 20 degree 
array<double>RandAngleJitter[RangDevision];
RandAngleJitter.fill(1,RangDevision,-10.0,-fillArrayJitter);
RandAngleJitter.fill(1+RangDevision/2,RangDevision,10.0,fillArrayJitter);

#term.print_line(fillArrayJitter);
loop int pp=1 until pp>RandAngleJitter.count()   # number of Trails
begin;
	#term.print_line(RandAngleJitter[pp]);
	RandAngleJitter[pp]=( RandAngleJitter[pp] + 90.0 )*pi_value/180.0;
	#term.print_line(RandAngleJitter[pp]);
	pp=pp+1;
end;

#term.print_line("Jitter  =  " +string(RandAngleJitter.count()));

#assiging random angles for all trials
array<double>StimulusType[0][0];
array<int>StimulusTypeCounter[Coherence.count()*Orient.count()];
StimulusTypeCounter.fill(1,0,0,0);
loop int pp=1 until pp>(Coherence.count()*Orient.count())  # number of Trails
begin;
	RandAngleJitter.shuffle();
	StimulusType.add(RandAngleJitter);
	pp=pp+1;
	
end;
#term.print_line("Type  =  "+string(StimulusType[1].count()));
# Initial graphic part for each frame
array<graphic_surface> trial_pictures[trial_duration_frames];
array<graphic_surface> AllTrial_pictures[num_trials*num_blocks][trial_duration_frames];
graphic_surface surface;
vsg::shape_generator shape = new vsg::shape_generator();  # line and shape generator
vsg::normalized_graphic norm_graph;
	# mask generator
vsg::gaussian_generator GaussianMask = new vsg::gaussian_generator( 0.0, 0.5*scrn_width);
vsg::normalized_graphic  mask = GaussianMask.generate( 2.0*scrn_width, 2.0*scrn_width );

# variable for circle shape creation 
array<double>thetaCircle[CircleDivision];
double fillArrayCircle=  pi_value/double(CircleDivision);# choosing semi random Theta angle.
thetaCircle.fill(1,0,0.0,fillArrayCircle); # add jittering for the each orientation
array<double>RandLocations_Final[trial_duration_frames][num_dots*CircleDivision][4];  # location of mini lines


# out put result file
output_file results_final = new output_file,
				results_final2 = new output_file;								
string file_name_output= logfile.subject() + "_results"+ ".txt",
		file_name_output2= logfile.subject() + "_MotionEnergy"+ ".txt";


#initialization of coherent and orientation
array<double>array_tmp[2];
loop int o=1 until o >Orient.count()   # number of orientations (here is two)
begin;
	loop int c=1 until c >Coherence.count()   #Coherence value for each trial
	begin;
		array_tmp[1]=double(c);#Coherence[c];
		array_tmp[2]=double(Orient[o]);
		Orient_Coherence.add(array_tmp);	
		int t=Orient_Coherence.count();
		#term.print_line( string(Orient_Coherence[t][2])+ " " +string(Orient_Coherence[t][1]));		
		c=c+1;
	end;
	o=o+1;
end;
	
# mirroring the naming if the fMRI_exp  =1
double mirror_value;
if fMRI_exp==1 then
	mirror_value=-1.0;
else
	mirror_value=1.0;
end;
	

################### Subroutine of Check_Location for MotionEnergy ########################################
sub int Check_AnglePos( double Xpos ,double Ypos)
begin
	int m = 0;
	if   Xpos>=0.0 && Ypos>=0.0 then	
		m=1;
	elseif Xpos<0.0 && Ypos>=0.0 then
		m=2;
	elseif Xpos<0.0 && Ypos<0.0 then
		m=3;
	elseif Xpos>=0.0 && Ypos<0.0 then
		m=4;
	end;
	return m;
end;



################### Subroutine of RandomPositions ####################
  # Generate non-overlapping starting locations
sub array<double, 2> RandomDots_position (int number_of_dots)
begin
	
   array<double>rand_locs[number_of_dots][2];
   array<int>num_dots_sector[0];
   
   double RDive=Rmax/double(KDivision),RRand,ThetaRand;
   int n=int(floor(double(number_of_dots/(KDivision*KDivision))));
	int tot_dot=0;
   # creating the sectors
	loop int j = 1; until j > KDivision
	begin
		num_dots_sector.add((2*j-1)*n);	
		tot_dot=tot_dot+num_dots_sector[j];
		j=j+1;
   end;
	#term.print_line(x_tmp);

   num_dots_sector[KDivision]=num_dots_sector[KDivision]+number_of_dots-tot_dot;	# the reaminging dots
   
	# we devide the square to 4 quadrants and create all coordinates in Radius and Theta. 
   #To do that, a quadrant would devide to equal surfaces which is (R^2*Pi)/(4*K^2). R radius and K the number of devisions of R.
	# the following will produce the XYs in first quadrant and add the adjust for other quadrants
	int nd=1;
	loop int j = 1; until j > KDivision
	begin
		loop int i = 1; until i > num_dots_sector[j]
		begin
			RRand=random()*RDive+double(j-1)*RDive; #choosing random radius
			ThetaRand= pi_value*double(random(1,360))/180.0; #choosing random radiuspi_value
			rand_locs[nd][1] =RRand *cos(ThetaRand);
			rand_locs[nd][2] =RRand *sin(ThetaRand);
			#term.print_line("Move  " +string(rand_locs[nd][1]) + "  " +string(rand_locs[nd][2]));
			i=i+1;
			nd=nd+1;
		end;	
		j = j + 1;
	end;
	#term.print_line("Done Here  ");
	rand_locs.shuffle();
	return rand_locs; # Return the 2-D array of X Y random values.
end;


################### Subroutine of Angle rotation ####################
sub double UpdateRandomDots_AngleOrient(double CohDotLables,double RandomAngle_jitter, double CoherOrient,double RandDirect)
begin
	double ThetaCalc,delta_theta_jitter;
	if CohDotLables==0.0 then
		delta_theta_jitter=0.0;#double(random(-1,1)*random(0,2))*pi_value/180.0;
		ThetaCalc=arccos(CoherOrient)+ RandDirect+delta_theta_jitter; #  random coherent 
		CurrentAngle=ThetaCalc;
	else	
		ThetaCalc=CohDotLables+ RandomAngle_jitter;
	end;
	return ThetaCalc;
end;

################### Subroutine of UpdatePositionCheck  &&  Check_Location ##############################
sub array<double, 2> UpdateAccording_AngleLife( array< double, 2 >& LocToUpdate,int NumCoherentDots,double OrientStim, 
		double RandDirect , int frame_number,array< double, 3 >& RandLocFrame)
		
begin
	array<double>Moved_Dots[0][0];
	array<double>QuadrantMotionEnergy[8];
	array<int> AliveDots[8];
	
	AliveDots.fill(1,0,0,0);
	QuadrantMotionEnergy.fill(1,0,0.0,0.0);
	
	## assigning the angle and orientation
	array<double>CohDotLables[num_dots];
   array<double>rand_theta_tmp[num_dots];
	double fillArray;
	double AngleUpdate,x_change,y_change,RCheck;

	Moved_Dots.assign(RandLocFrame[random(1,1000)]); #assignng a random position for all dots for updating if needed 
	
	CohDotLables.fill(1,0,0.0,0.0); # all angle are zeros
	
	double fillShift = (pi_value/180.0)*double(random(0,60));  # generating a random shift for all non coherent orientations 
	
	if num_dots-NumCoherentDots>0 then
		fillArray=2.0*pi_value/double(num_dots-NumCoherentDots);# the fillarray for adding to the array
		CohDotLables.fill(1,num_dots-NumCoherentDots,0.0,fillArray); # none coherent angles assigned to the dots 	
	end;
	
	CohDotLables.shuffle(); # shuffling all dots  (the coherent dots has 0.0 value for angle
	CohDotLables.shuffle();
	int DeadCohDot=0;
	
	loop int m = 1 until m > num_dots
	begin
		
		# generating the usable angleupdate
		AngleUpdate = UpdateRandomDots_AngleOrient(CohDotLables[m],fillShift,OrientStim,RandDirect);
		x_change=move_amount * cos(AngleUpdate);
		y_change=move_amount * sin(AngleUpdate)  ;
		
		# Updating the Location of dots
		LocToUpdate[m][1]= LocToUpdate[m][1] +   x_change  ;
		LocToUpdate[m][2]= LocToUpdate[m][2] +   y_change ;  
		
		# cheking the life time
		lifetime_matrix[m]=lifetime_matrix[m]-1;
		# cheking the location
		RCheck=sqrt(LocToUpdate[m][1]*LocToUpdate[m][1] +LocToUpdate[m][2]*LocToUpdate[m][2]);
		int SectorNum=0;

		# check the life time and the location
		if RCheck >Rmax ||  lifetime_matrix[m]==0  then
			lifetime_matrix[m]=lifetime_matrix_fixed[m];
			if CohDotLables[m]==0.0 then;
				DeadCohDot=DeadCohDot+1;
			end;
		else
			
			if RCheck>Rmax/2.0 then
				SectorNum=Check_AnglePos(LocToUpdate[m][1],LocToUpdate[m][2])+4;
			else
				SectorNum=Check_AnglePos(LocToUpdate[m][1],LocToUpdate[m][2]);
			end;

			AliveDots[SectorNum]=AliveDots[SectorNum]+1;
			QuadrantMotionEnergy[SectorNum]=QuadrantMotionEnergy[SectorNum]+ x_change;	
			Moved_Dots[m].assign(LocToUpdate[m]);
			
		end;	

		# RandLocat_CircleGeneration
		loop int j = 1 until j > CircleDivision
		begin
			int dot_num= (m-1)*CircleDivision + j ;
			RandLocations_Final[frame_number][dot_num][1]  =  Moved_Dots[m][1] +  dot_radius * cos(thetaCircle[j]+pi_value);
			RandLocations_Final[frame_number][dot_num][2]  =  Moved_Dots[m][2] +  dot_radius * sin(thetaCircle[j]+pi_value);
			RandLocations_Final[frame_number][dot_num][3]  =  Moved_Dots[m][1] +  dot_radius * cos(thetaCircle[j]);
			RandLocations_Final[frame_number][dot_num][4]  =  Moved_Dots[m][2] +  dot_radius * sin(thetaCircle[j]);
			#term.print_line("  x1  " +string(x1) + "  y1  " +string(y1) +"  x2  " +string(x2) + "  y2  " +string(y2) );
			j=j+1;
		end;
		
		
		m=m+1;
	end;
	
	
	#arithmetic_mean
	loop int ss = 1 until ss > 8
	begin
		if AliveDots[ss]>0 then
			TrailMotionEnergy_frame[frame_number][ss]=QuadrantMotionEnergy[ss]/double(AliveDots[ss]);
		end;
		ss=ss+1;
	end;
	
	return Moved_Dots;
end;


################### Subroutine of UpdateRandomPositions ####################
sub array<graphic_surface,1>  UpdateRandomDots_position(int trial_type,array< double, 3 >& RandLocFrame)
begin

	# initial array for each trial
	array<int>coherent_tmp[0][0];
	double RandomDirection=0.0;
	double AllColor= 100.5;
	
	#rand_locations_frame[1]=RandomDots_position(num_dots); # initialzed location for dots
	rand_locations_frame[1]=RandLocFrame[random(1,1000)]; # initialzed location for dots
	# assiginging a life time to each dot
	# it has a life counter on the first coulumn and 
   loop int m = 1 until m > num_dots
	begin
		lifetime_matrix[m]=random(1,MaxLifeTime);
		m=m+1;
	end;
	lifetime_matrix_fixed.assign(lifetime_matrix);
	
	int pp=0;
	double Coh=Orient_Coherence[trial_type][1];
	
	double OrientStim=Orient_Coherence[trial_type][2];  # the orientation
	#term.print_line( Coh*OrientStim ); 
	double CoherentStim= Coherence[int(Coh)]/100.0; # the coherent
	int NumCoherentDots= int(floor(double(num_dots)*CoherentStim));# number of coherent dots 
	
	CurrentCoherent=Coherence[int(Coh)]*OrientStim;
	
	if CurrentCoherent <0.0 then 
		pp=0;
	else
		pp=Coherence.count();
	end;
	
	pp=pp+int(Coh);
	
	if CheckFirstTrial==0 then
		RandomDirection=StimulusType[1][1];
	else
		StimulusTypeCounter[pp]=StimulusTypeCounter[pp]+1;
		#term.print_line("PP Num : "+string(StimulusTypeCounter[pp])+ "   "+ string(pp));
		RandomDirection=StimulusType[pp][StimulusTypeCounter[pp]];
		
	end;
	
	rand_locations_frame[1]=UpdateAccording_AngleLife(rand_locations_frame[1], NumCoherentDots,OrientStim,RandomDirection ,1, RandLocFrame);	
	#int StartTime1=clock.time();
	loop int f = 2 until f > trial_duration_frames
	begin
		# controling if any dots has changed
		rand_locations_frame[f]=UpdateAccording_AngleLife(rand_locations_frame[f-1], NumCoherentDots,OrientStim,RandomDirection ,f, RandLocFrame);
		
		shape.add_line_set( RandLocations_Final[f-1],  line_width, 1.0 , 1.0, 0.0 ); # generating small lines
		norm_graph = shape.generate( 2.0*scrn_width, 2.0*scrn_width ); # generating norm_graph from shapes
		norm_graph.mask2(mask );  # mask the shape with guassian mask
		norm_graph.multiply_add(2.0,-1.0); 
		surface=norm_graph.create_graphic(AllColor, AllColor, AllColor, AllColor, AllColor, AllColor );
		#surface.flip_horizontally();
		trial_pictures[f-1]=surface;	
		shape.clear( );
		f=f+1;
	end;
	#term.print_line("TimeCalc : "+string(clock.time()-StartTime1));
	return trial_pictures;
end;


################### Subroutine of RandomDots ##########################

sub LOGFILE_Writing(array< double, 1 >& log_response_stmuli)
begin	
	if  log_response_stmuli[1]==1.0 then
		results_final.open ( file_name_output, true);
		results_final.print("################# Begining ########## \n\n" );
		results_final.print("##  L/R   Left/Right \n" );
		results_final.print("##  Coh   Coherent percentage \n" );
		results_final.print("##  A   Angle  \n" );
		results_final.print("Trial\tCoh\tL/R\tRT\tC/W\tA\tSt\tRt  \n" );
		results_final.print("  \n" );
	end;
	string A="";
	loop  int k=1 until k> log_response_stmuli.count() 
	begin
		A=A+string(log_response_stmuli[k]) + "\t"  ;
		k=k+1;
	end;
	results_final.print( A + "\n" );
end;

################### Subroutine of Final Performance ##########################
sub LOGFILE_ResultWriting
begin
	array<double>ResultPerformance[2*Coherence.count()][5];
	
	loop  int k=1 until k> TotalLogData.count() 
	begin
		loop  int p=1 until p> Coherence.count() 
		begin
			if TotalLogData[k][2]>0.0 && (abs(TotalLogData[k][2])==Coherence[p]) then
				ResultPerformance[p][1]=Coherence[p];
				if TotalLogData[k][5]>0.0 then # if the trial is correct
					ResultPerformance[p][2]=ResultPerformance[p][2]+TotalLogData[k][5];  # summing correct responses
					ResultPerformance[p][3]=ResultPerformance[p][3]+TotalLogData[k][4];  # summing correct RT
					ResultPerformance[p][5]=ResultPerformance[p][5]+1.0;
				end;
				ResultPerformance[p][4]=ResultPerformance[p][4]+double(1);
			elseif TotalLogData[k][2]<0.0 && (abs(TotalLogData[k][2])==Coherence[p]) then
				int q=p+Coherence.count();
				#term.print_line(q);
				ResultPerformance[q][1]=-1.0*Coherence[p];
				if TotalLogData[k][5]>0.0 then # if the trial is correct
					ResultPerformance[q][2]=ResultPerformance[q][2]+TotalLogData[k][5]; # summing correct responses
					ResultPerformance[q][3]=ResultPerformance[q][3]+TotalLogData[k][4]; # summing correct RT
					ResultPerformance[q][5]=ResultPerformance[q][5]+1.0;
				end;
				ResultPerformance[q][4]=ResultPerformance[q][4]+1.0;
			end;
			p=p+1;
		end;	
		k=k+1;
	end;
	loop  int p=1 until p> ResultPerformance.count() 
	begin	
		if ResultPerformance[p][2]>double(0) then
			ResultPerformance[p][3]=ResultPerformance[p][3]/ResultPerformance[p][2];
		else
			ResultPerformance[p][3]=ResultPerformance[p][2];
		end;
		ResultPerformance[p][2]=ResultPerformance[p][2]/ResultPerformance[p][4];
		p=p+1;
	end;
	#term.print_line(log_response_stmuli.count());

	results_final.print("  \n" );
	results_final.print("################# Result Part ########## \n\n" );
	results_final.print("Cond\tPerf\tRT\tTotal\tCor\n" );
	results_final.print("  \n" );
	string A="";
	loop  int p=1 until p> ResultPerformance.count() 
	begin
			A="";
			loop  int k=1 until k> ResultPerformance[1].count() 
			begin
				A=A+string(ResultPerformance[p][k]) + "\t"  ;
				k=k+1;
			end;
			results_final.print( A + "\n" );
		p=p+1;
	end;
	
	# adding the motion enery data to the Text file
	results_final2.open ( file_name_output2, true);
	loop  int t=1 until t>TrailMotionEnergy_trails.count()
		begin
			loop  int ss=1 until ss> 8
			begin
				A="";
				loop  int f=1 until f> trial_duration_frames
				begin
					A=A+string(double(int(TrailMotionEnergy_trails[t][f][ss]*10000.0))/10000.0) + ","  ;
					f=f+1;
				end;
				ss=ss+1;
				results_final2.print( A + "\n" );
			end;
		t=t+1;
	end;
end;


#################fixation subrutine reset ############################
sub FixationDot_rest
begin
	EllipseShape.set_dimensions( FixEllipse, FixEllipse );
	EllipseShapeGreen.set_dimensions( FixEllipse, FixEllipse );
	EllipseShapeGreen.set_color(55, 255, 0, 200 );
	EllipseShapeGreen.redraw();
	EllipseShape.set_color( 255, 0, 0, 200 );
	EllipseShape.redraw();
#	pic.insert_part(2,EllipseShape,0, 0);
end;

#=====================================
sub WAIT( int i )
begin
	int start = pulse_manager.main_pulse_count();
   loop
   until pulse_manager.main_pulse_count() >= start + i 
   begin    
   end; 
end;

###################### Generating InterTrial Intervals and Dummy Trials ###########################

sub  InterTrial_JitterGen
begin
	int ii=0;
	string x_tmp="";
	loop int pp=1 until pp>MSeq.count()   # number of Trails
	begin;
		if pp<3*MSeq.count()/4 then
			ii=0;
		elseif pp>=3*MSeq.count()/4   then  #&& pp<=(5*num_trials)/6
			ii=1;
		#else
		#	ii=2;
		end;
		#ii=0;
		InterTrialDealy[pp]=ii;
		 x_tmp=x_tmp+ " " +string(InterTrialDealy[pp]);	
		pp=pp+1;
	end;
	#term.print_line( "Intervals  = " + x_tmp);
end;
################### Subroutine of RandomDots ##########################
sub array<double, 1> RandomDots_present( array<graphic_surface,1> &FramePic,
			int previous_counting, int trial_type, int trial_number,int real_trial_number)
begin
	array<double>log_data[8];
	int button_diff = 0,  
		NumRespToDo=0,
		num_button_diff = 0,  
		ReactionTime=0,
		button = 0,
		accepted_duration=accepted_trial_duration,
		i_tmp=0,
		lonResp=0;
	string LogTime="";		
	
	# number of button presses
	if MSeqPerm[trial_number] >6 then
		NumRespToDo=8;
		lonResp=2;
	else
		NumRespToDo=3;
		lonResp=1;
	end;
	
	log_data[1]=double(real_trial_number); # Log Counter
	
	double tmp_cod=0.0,tmp_mir=0.0,Stim_main_code=0.0;
	string direction_tmp="";
	if CurrentCoherent*mirror_value>0.0 then
		direction_tmp="up";
		tmp_cod=20.0
	else
		direction_tmp="down";
		tmp_cod=10.0;
	end;
	   # registering stimuli type
	double stim_type=mirror_value*Orient_Coherence[trial_type][2] * Coherence[int(Orient_Coherence[trial_type][1])];
	log_data[2]=stim_type;
	
	

	double port_code_var=Orient_Coherence[trial_type][1]*10.0+ tmp_cod*10.0+double(lonResp);  # building the code output

	Stim_main_code=floor(CurrentAngle*180.0/pi_value);
	
	#term.print_line(string(NumRespToDo) +"   " +string(stim_type) +"   " + string(MSeqPerm[trial_number]));
	
   main_event.set_event_code( string(Stim_main_code));
	log_data[6]=Stim_main_code;
	#DummyEvent.set_event_code( string(stim_type));
	#term.print_line(string(Stim_main_code)+"   " +string(floor(CurrentAngle*180.0/pi_value)) );
	#main_trial.present();  # adding to the log file
	# Generating the random dots
	int reaction_check=0;
	int stop_trial=0,StimTime=0,RespTime=0;
	int i = 1,
		StartFrame=0;
	int write_stop_trial=0;
	int passed_duration_time=0;
	#term.print_line(pic.part_count());
	DummyEvent.set_event_code( string(port_code_var));# writing the green fixation
	DummyTrial.present();
	pic.insert_part(1,FramePic[i],0, 0);
	main_event.set_port_code(port_code_var);
	#extrae_vent.set_event_code(string(port_code_var));
	main_trial.present();
	LastStimFirst=stimulus_manager.last_stimulus_data();
	StimTime=LastStimFirst.time()-ExpStartVolTime;
	LastStim = stimulus_manager.last_stimulus_data();
	i_tmp=i;
	DummyEvent.set_event_code( string(999));# writing the green fixation
	
   loop 
      int StartT = clock.time();
   until 
       passed_duration_time>accepted_duration || button_diff>=NumRespToDo  # Show picture for 500 msec
   begin

		if i > i_tmp  then
			pic.remove_part(1);
			pic.insert_part(1,FramePic[i],0, 0);
			pic.present(); # here showing the pic
			StartFrame = clock.time();
			i_tmp=i;
		end;

		button_diff = response_manager.response_data_count()- previous_counting; # assessing last response (1 or 2)
		if button_diff==1 && write_stop_trial==0 then
			LastResp = response_manager.last_response_data();

			pic.remove_part(2);
			write_stop_trial=write_stop_trial+1;
			pic.insert_part(2,EllipseShapeGreen,0, 0);
			#DummyTrial.present();
			pic.present(); # 
			StimEvent.set_event_code( string(999));  # Stimulus Disapeared timing
			TrialName.present();
			
			button=response_manager.last_response();

			ReactionTime=LastResp.time()-LastStim.time( );
			RespTime=LastResp.time()-ExpStartVolTime;
			
		end;


		passed_duration_time=clock.time()-StartT;
		if passed_duration_time<1800 && button_diff>0 && reaction_check==0 then
			if NumRespToDo >3 then   # if they pressed one button it adds up the timing till 2.2 sec
				accepted_duration=passed_duration_time+2200; # adding more time 
			else
				accepted_duration=passed_duration_time+1000;
			end;
			reaction_check=reaction_check+1;
		end;
		
		
		int InterFrameDelay=clock.time()-StartFrame;
		if  InterFrameDelay >= wait_interval_time_btwenn_frame then
			i = i + 1;
			pic.remove_part(1);
			pic.insert_part(1,EllipseShapeBlack,0,0);
			if i> trial_duration_frames-2 then
				i=1;i_tmp=0;
			end;
		end;
		
		#term.print_line("Frame= " +string(i) + "  Delay = "+string (InterFrameDelay));
		LogTime=LogTime+"\n"+"Frame= " +string(i) + "  Delay = "+string (InterFrameDelay) + "  CurrDur= " 
				+string(passed_duration_time)+ "  ChagDur= " +string(accepted_duration)+"  Dur= " +string(accepted_trial_duration) ;
   end;	
	pic.remove_part(2);
	pic.remove_part(1);
	pic.insert_part(1,EllipseShape,0, 0);
	pic.present();

	# if there is no button pressed
	if button_diff==0 then
		DummyEvent.set_event_code( string(999));
		DummyTrial.present();
	end;
	StimEvent.set_event_code( string(888));  # Stimulus Disapeared timing
	TrialName.present();
	
	#term.print_line("Stm= " +string(stim_type) );
	int CorrectResponse=0;
   if ReactionTime!=0  then
		if stim_type<0.0 && button ==up then
			 CorrectResponse=1;
		elseif stim_type>0.0 && button ==down then
			 CorrectResponse=1;
		end;
	end;
	#term.print_line(LogTime);
	log_data[3]=double(button);
	log_data[4]=double(ReactionTime);
	log_data[5]=double(CorrectResponse);
	log_data[7]=double(StimTime)/1000.0;
	log_data[8]=double(RespTime)/1000.0;
	return log_data
end;


###################### Main Programm ###########################

# generating the whole matrix of the experiment         
 #randomizing the X and Y 

 FixationDot_rest();
array<graphic_surface>FramePic_trial[0];
array<double>RandLovAllFrames[0][0][0];
array<double>TrialLogData[0];
int accurate_response=0;
int missed_response=0;
int wait_time,volnumber;

# Generationg all random locations for all the frames
set_pcl_exclusive_mode(true);
loop int f = 1 until f > initial_random_frames
begin
	RandLovAllFrames.add(RandomDots_position(num_dots));
	f=f+1;		
end;
# cheking the speed for calculation
int TrialCreationTiming=clock.time(),timecalc;
; #updating rand_locations_frame and RandLocations_Final array
FramePic_trial.assign(UpdateRandomDots_position(1,RandLovAllFrames));  # Assign dots to graphic object
set_pcl_exclusive_mode(false);
timecalc= clock.time()-TrialCreationTiming;  # the time to calculate one trial
#term.print_line("Time of a trial= " +string(timecalc) );
CheckFirstTrial=1;

int NumSmallJitter=VolumeTR/(MSeq.count()*num_blocks);
array<int>TrialSmallJitter[MSeq.count()*num_blocks];

TrialSmallJitter.fill(1,0,0,NumSmallJitter);
TrialSmallJitter.shuffle();

InterTrial_JitterGen();

int resp_count=0;

DummyEvent.set_event_code( "00000");
DummyTrial.present();
term.print_line("EXPTimeStart= " +string(clock.time()-TrialCreationTiming) );

#term.print_line(xtmp);

timing.set_caption("Get Ready" );
timing.redraw();
if mirror_value<0.0 then
	timing.flip_horizontally();
end;
long_PAUSE.present();
#First Picture 


WAIT(1); 
ExpStartVolTime=pulse_manager.main_pulse_time(1);

pic.present();
# Count down padding volumes
WAIT(padding_volumes);


InterTrialDealy.shuffle();
Orient_Coherence.shuffle();

pic.present();
int real_trial_counter=0;
int t=0;

loop int b=1 until b > num_blocks# Orient_Coherence.count();   # number of Trails
begin;
	int start_block=clock.time();
	#InterTrialDealy.shuffle();
	int DumCounter=1, InterTrial=0;
	#Orient_Coherence.shuffle();
	int Block_trial_counter=0;
	loop int tt=1 until tt > num_trials_perBlock # Orient_Coherence.count();   # number of Trails
	begin;
		timecalc= clock.time();
		t=t+1;
		if MSeq[t]>0 then
			real_trial_counter=real_trial_counter+1;
			Block_trial_counter=Block_trial_counter+1;
			set_pcl_exclusive_mode(true);
				#updating rand_locations_frame and RandLocations_Final array
			#FramePic_trial.assign(UpdateRandomDots_position(MSeq[t],RandLovAllFrames));  # Assign dots to graphic object
			FramePic_trial.assign(UpdateRandomDots_position(MSeq[t],RandLovAllFrames));  # Assign dots to graphic object
			set_pcl_exclusive_mode(false);

			WAIT(1);
			#wait_time=VolumeTR +TrialSmallJitter[(b-1)*MSeq.count() + t] + random(1,50)-timecalc+InterTrialDealy[t]*VolumeTR;
			wait_time=VolumeTR +TrialSmallJitter[t] + random(1,50)+InterTrialDealy[t]*VolumeTR ;# -timecalc
			
			DummyEvent.set_event_code("0000"+ string(wait_time));
			DummyTrial.present();
			wait_interval(wait_time);  # inter trial interval

			resp_count =response_manager.response_data_count();
			volnumber=pulse_manager.main_pulse_count();
			
			InterTrial=clock.time()-timecalc;  # tim to calculate 
			if InterTrial<4400 then # minimum intertrial interval
				wait_interval(abs(4400-wait_time));
				term.print_line(InterTrial);
			end;
			
			#TrialLogData.assign(RandomDots_present(FramePic_trial,resp_count,MSeq[t],real_trial_counter)); # feedbavck from presenting
			TrialLogData.assign(RandomDots_present(FramePic_trial,resp_count,MSeq[t],t,real_trial_counter)); # feedbavck from presenting
			
			TotalLogData[real_trial_counter]=TrialLogData;
			#term.print_line(wait_time);
			#term.print_line(string(TrialLogData[7]*1000.0-double(pulse_manager.main_pulse_time(volnumber))));
			
			#pic.insert_part(1,EllipseShape,0, 0);
			#pic.present();
				
			LOGFILE_Writing(TrialLogData);
			TrailMotionEnergy_trails[ real_trial_counter ]= TrailMotionEnergy_frame; #updating motion energy matrix
			
			accurate_response=accurate_response+int(TrialLogData[5]); #counting the correct responses
			if TrialLogData[3]==0.0 then
				missed_response=missed_response+1; #counting the missed responses
			end;
			RandLovAllFrames.shuffle();
		else	
			wait_time=VolumeTR +TrialSmallJitter[ t] + random(1,50)+InterTrialDealy[t]*VolumeTR;
			wait_interval(wait_time);  
			main_event.set_event_code( string(000));
			main_event.set_port_code( 100);
			main_trial.present();
			if DumCounter<num_trials_percond then
				DumCounter=DumCounter+1;
			end;
			wait_interval(trial_duration);
			wait_interval(wait_time+random(10,100));
		end;
		tt=tt+1;	# trial per block 
	end;
	# Block information Check
	wait_time=VolumeTR-random(300,500);
	wait_interval(1000+wait_time);  # time before feedback
	
	int BlockPercentage=(b*100)/num_blocks;
	
	int AccuracyData=int(floor(double(100*accurate_response/Block_trial_counter)));
	timing.set_caption(" <b>Accuracy = " +string(AccuracyData) +"%</b>" +" \n\n <b>Missed Responses = "  
			+ string(missed_response) +"</b>");
	timing1.set_caption(" Percent completed  = " +string(BlockPercentage) +"%" );
	timing1.redraw();
	if mirror_value<0.0 then
		timing1.flip_horizontally();
	end;
	timing.set_formatted_text( true );
	timing.redraw();
	if mirror_value<0.0 then
		timing.flip_horizontally();
	end;
	long_PAUSE.present();
	EllipseShape.set_color( 155, 155, 155, 200 );
	EllipseShape.redraw();
	pic.present();
	
	wait_interval(wait_time+7000);  # time after feedback
	
	EllipseShape.set_color( 255, 0, 0, 200 );
	EllipseShape.redraw();
	#pic.insert_part(2,EllipseShape,0, 0);
	pic.present();
	accurate_response=0;
	missed_response=0;
	b=b+1;
	start_block=clock.time()-start_block;
	#term.print_line("BlockTime= " +string(start_block) );
end;	
WAIT(padding_volumes);# change it to 10
LOGFILE_ResultWriting();  #writing the final log 
term.print_line("EXPTimeEnd= " +string(clock.time()-TrialCreationTiming) );


