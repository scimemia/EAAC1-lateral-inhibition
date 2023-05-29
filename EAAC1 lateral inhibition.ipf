Proc AverageSweeps(wn,chan,rnum,sweepvect)
	String wn="AS052923A",sweepvect="sweeps"
	Variable rnum=1,chan=1,disp=0
	Prompt wn, "File name (no path, no .ibw):"
	Prompt chan, "Channel number: "
	Prompt rnum, "Run number: "
	Prompt sweepvect, "Wave containing sweep numbers:"
	Silent 1
	variable numtimes=numpnts($sweepvect),jump,i,blip,mult,trodes
	string dest
	jump = $wn[2*rnum-1]*256												//Set pointer in file
	trodes=$wn[jump+6]													//Number of channels
	if(trodes<=0)															//Backwards compatibility
		trodes = 1
	endif
	dest=wn+"_ch"+num2str(chan)+"_"+num2str(rnum)+"_ave_"+num2str($sweepvect[0])+"_"+num2str($sweepvect[numpnts($sweepvect)-1])
	Make /O/N=(256*$wn[jump+8]-3) $dest; $dest=0							//[jump+8] is blocks per run
	SetScale/P x 0,($wn[jump+5]/100),"", $dest
	i=0
	do
		$dest += $wn[p+jump+256+(trodes*($sweepvect[i]-1)+chan-1)*256*$wn[jump+8]]
		i+=1
	while(i<numtimes))
	$dest/=numtimes
	Wavestats/Q/R=(0,2) $dest
	DoAlert 1, "Subtract baseline?"
		if (V_flag==1)
			$dest-=V_avg							//subtract baseline
		endif							
	//if(AV_Alert)										//Allows programmer to turn off feedback for use in macros
		DoAlert 1, "Make a new graph?"
		if (V_flag==1)
			display $dest;ModifyGraph lsize=0.5
			ShowInfo
		else
			DoAlert 1, "Append wave to existing graph?"
			if (V_flag==1)
				AppendtoGraph $dest
			endif
		endif
		print dest
	//endif
End Macro

//AnalyzeEPSC_AS uses different ways to calculate the decaying phase of EPSCs: mono- and bi-exponential fits, weighted time constants, time to 50 and 10%, 
//tcourse over a 600 ms time window. It works with different initial guesses for NMDA and AMPA EPSCs.

Proc AnalyzeEPSC_AS(wn,int,stlat,io, intw, pn)
	String wn="AS000005A_2_ave_1_10_norm"
	Variable int=0.2,stlat=100, io=1, intw=600, pn=1
	Prompt wn, "Wave name",popup WaveList("*",";","WIN:")
	Prompt int,"Sampling rate"
	Prompt stlat, "Stimulus latency"
	Prompt io, "Inward (-1) or Outward (1)"
	Prompt intw, "Integration window (ms)"
	Prompt pn, "Pulse number"
	Silent 1
	
	string singfit=wn+"_sx_p"+num2str(pn),dblfit=wn+"_dx_p"+num2str(pn)
	variable peakloc,peakval,peaklat,baseval,startfit,halfdecaytime,tendecaytime,tcourse,frac1,frac2,rt02,rt08,rt0208,t50f,t50b,errort50,t45,t55,HW1,HW2,HW,tau
	Duplicate/O $wn tvsm
	Setscale/P x 0, int,"", tvsm
	WaveStats /Q/R=(stlat-20,stlat-10) tvsm; baseval=V_avg
	Smooth 30, tvsm
	tvsm*=io																				
	FindPeak/B=10/P/Q/R=(xcsr(A),xcsr(B)) tvsm; 
	peakloc=round(V_PeakLoc)*int; peakval=tvsm[round(V_PeakLoc)]
	EdgeStats/B=10/L=(peakval,baseval)/Q/R=(peakloc,xcsr(B))/T=30 tvsm							//Finds the 10-90% of the decaying phase you want to fit
	startfit = V_EdgeLoc1
	peaklat=peakloc-stlat
	halfdecaytime=V_EdgeLoc2-peakloc
	tendecaytime=V_EdgeLoc3-peakloc
	//Duplicate/O $wn tvsm$singfit,$dblfit;$singfit=NaN;$dblfit=NaN
	tvsm-=baseval																			//subtract baseline
	tvsm/=abs(peakval)																		//normalizes the wave by the peak
	Make/O/N=3 sing_coef,sing_sigma;Make/O/N=5 dbl_coef,dbl_sigma
	//K0=0;K1=2*io;K2=.01																		//initial guesses for exp fit + constrain to zero baseline
	//CurveFit/Q/H="100" exp tvsm(startfit,xcsr(B)) /D=$singfit
	//Make/O/N=4 W_coef,W_sigma
	//sing_coef=W_coef;sing_sigma=W_sigma;
	//tau=1/K2
	//if (io==1)
	//	K0=0;K1=0.5*sing_coef[1];K2=0.7*sing_coef[2];K3=0.7*sing_coef[1];K4=1.2*sing_coef[2]		//initial guesses for dblexp fit on outward NMDA-like EPSCs
	//endif
	//if (io==-1)																				//initial guesses for dblexp fit on inward AMPA-like EPSCs
	//	K0=0;K1=0.5*sing_coef[1];K2=0.7*sing_coef[2];K3=0.7*sing_coef[1];K4=1.2*sing_coef[2]		//multiplying factors for AMPA EPSCs: 12;1;7;0,85
	//endif
	//CurveFit/Q/G/L=5000/H="10000" dblexp tvsm(startfit,xcsr(B)) /D=$dblfit
	//dbl_coef=W_coef;dbl_sigma=W_sigma;
	//frac1=abs(dbl_coef[1])/(abs(dbl_coef[1])+abs(dbl_coef[3]));frac2=abs(dbl_coef[3])/(abs(dbl_coef[1])+abs(dbl_coef[3]))
	Wavestats/Q/R=(V_PeakLoc*int,V_PeakLoc*int+intw) tvsm
	tcourse=abs(V_avg*intw)																	//there is no need to divide by peakval because tvsm has already been normalized by pkval
	FindLevel/B=3/Q/R=(xcsr(A),peakloc) tvsm, 0.2
	rt02=V_LevelX; //print V_levelX
	FindLevel/B=3/Q/R=(xcsr(A),peakloc) tvsm, 0.8
	rt08=V_LevelX; //print V_levelX
	rt0208=rt08-rt02																			//calculates the 20-80% rise time
	FindLevel/B=5/Q/R=(stlat+0,peakloc) tvsm, 0.5
	HW1=V_LevelX
	FindLevel/B=5/Q/R=(peakloc,xcsr(B)) tvsm, 0.5											//calculates the half width
	HW2=V_LevelX
	HW=HW2-HW1
	FindLevel/B=5/Q/R=(peakloc,xcsr(B)) tvsm, 0.55												//these 2 commands set the 50+/-5% time window where to look for the t50
	t55=V_LevelX
	FindLevel/B=5/Q/R=(xcsr(B), peakloc) tvsm, 0.45
	t45=V_LevelX
	FindLevel/B=1/Q/R=(t55,t45) tvsm, 0.5														//find the level in the left-to-right direction
	t50f=V_LevelX
	FindLevel/B=1/Q/R=(t45, t55) tvsm, 0.5														//find the level in the right-to-left direction
	t50b=V_LevelX
	errort50=t50b-t50f
	printf "Single Exp fit:\r\ttau =%5.1f ms.\r",1/sing_coef[2]
	printf "Double Exp fit:\r\ttaufast =%5.1f ms (%2.0f%%)\r\ttauslow =%5.1f ms (%2.0f%%)\r",1/dbl_coef[4],100*frac1,1/dbl_coef[2],100*frac2
	printf "Peak Amplitude = %g pA\r",peakval-baseval
	printf "Peak Latency = %g ms\r", peaklat
	printf "Rise Time (20-80) = %g ms\r", rt0208
	printf "Exp decay (tau) = %g ms\r",tau
	printf "Half decay time = %g ms\r",halfdecaytime
	printf "Half decay time error = %g ms\r", errort50
	printf "Ten decay time = %g ms\r",tendecaytime
	printf "tcourse = %g ms\r",tcourse
	printf "half width = %g ms\r",HW
	Killwaves sing_coef,sing_sigma,dbl_coef,dbl_sigma
	DoAlert 1, "Append fit waves?"
	if (V_flag==1)
		//SetScale/P x 0,int,"", $singfit,$dblfit
		//AppendToGraph tvsm,$singfit,$dblfit
		//ModifyGraph rgb(tvsm)=(43520,43520,43520);DelayUpdate
		//ModifyGraph rgb($singfit)=(0,12800,52224);DelayUpdate
		//ModifyGraph rgb($dblfit)=(0,52224,0)
	else
		DoAlert 1, "Kill fit waves?"
		if (V_flag==1)
			//Killwaves $singfit,$dblfit
		endif
	endif
	
EndMacro