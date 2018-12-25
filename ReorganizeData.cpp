/*
  ReorganizeData.cpp: This program imports the xU_t().csv files
  output by the main code and splits each file into multiple
  sub-files, one per field variable. The objective here
  is to minimize the instantaneous memory consumption of the
  accompanying Matlab plotting schemes.
*/

#include<stdio.h>
#include<math.h>
#include<vector>
#include<iostream>


//Mesh resolution: code reads it in from size_params.csv
int Mx;
int My;

//fluid parameters
const double rsh = 1.4;

typedef std::vector<std::vector<double> > D2_double;
typedef std::vector<std::vector<std::vector<std::vector<double> > > > D4_double;

D2_double Allocate_2D_double(int d1, int d2)
{
  std::vector<std::vector<double> > output;
  output.resize(d1);
  for (int A = 0; A < d1; A = A+1)
    {
      output[A].resize(d2);
    }
  return output;
}
D4_double Allocate_4D_double(int d1, int d2, int d3, int d4)
{
  std::vector<std::vector<std::vector<std::vector<double> > > > output;
  output.resize(d1);
  for (int A = 0; A < d1; A = A + 1)
    {
      output[A].resize(d2);
      for (int B = 0; B < d2; B = B + 1)
	{
	  output[A][B].resize(d3);
	  for (int C = 0; C < d3; C = C + 1)
	    {
	      output[A][B][C].resize(d4);
	    }
	}
    }
  return output;
}

void transcribe(int step)
{
  //Write the subfiles for one specific substep, presumably out of 9.
  std::string filename;
  char bufferStep [50];
  int nstep = sprintf(bufferStep,"%d",step);
  printf("Preparing to import code output, step %d\n", step);
  FILE*input;
  std::string FileName = "xU_t";
  FileName.append(bufferStep);
  FileName.append(".csv");
  input = fopen(FileName.c_str(),"r");
  //Import the entire set of information
  double messenger;
  double garbage;

  D2_double xrec = Allocate_2D_double(Mx,My);
  D2_double yrec = Allocate_2D_double(Mx,My);
  
  D2_double rho_local = Allocate_2D_double(Mx,My);
  D2_double vx_local = Allocate_2D_double(Mx,My);
  D2_double vy_local = Allocate_2D_double(Mx,My);
  D2_double p_local = Allocate_2D_double(Mx,My);
  double en_local;
  printf("Grabbing from xU_t%d.csv\n",step);
  for (int I = 0; I < Mx; I++)
    {
      for (int J = 0; J < My; J++)
	{
	  //get xrec and yrec entries
	  fscanf(input,"%lf,", &messenger);
	  xrec[I][J] = messenger;
	  fscanf(input,"%lf,", &messenger);
	  yrec[I][J] = messenger;
	  
	  //Now, get rho,u,v,specific energy
	  fscanf(input,"%lf,",&messenger);
	  rho_local[I][J] = messenger;
	  fscanf(input,"%lf,", &messenger);
	  vx_local[I][J] = messenger;
	  fscanf(input,"%lf,", &messenger);
	  vy_local[I][J] = messenger;
	  fscanf(input,"%lf,", &messenger);
	  en_local = messenger * rho_local[I][J];
	  p_local[I][J] = (rsh-1.0) * (en_local - 0.5*rho_local[I][J]*(vx_local[I][J]*vx_local[I][J] + vy_local[I][J]*vy_local[I][J]));

	  //skip velocity entry + pressure entry
	  fscanf(input,"%lf,", &messenger);
	  garbage = messenger;
	  fscanf(input,"%lf,", &messenger);
	  garbage = messenger;
	}
    }
  if (step == 0)
    {
      //If inspecting the first output file, grab the geometry. This only needs to happen once.
      FILE*fileGeo = fopen("XY.csv","w");
      for (int I = 0; I < Mx; I++)
	{
	  for (int J = 0; J < My; J++)
	    {
	      fprintf(fileGeo,"%8.7f, %8.7f\n",xrec[I][J],yrec[I][J]);
	    }
	}
      fclose(fileGeo);
    }
  printf("Preparing to transcribe the data to separate sub-files at step %d\n", step);
  //Transcribe the data to a set of output files particular to the current step.
  FILE*out_rho;
  FileName = "rho_t";
  FileName.append(bufferStep);
  FileName.append(".csv");
  out_rho = fopen(FileName.c_str(),"w");
  //printf("passed first Filename overwirte\n");
  FILE*out_vx;
  FileName = "vx_t";
  FileName.append(bufferStep);
  FileName.append(".csv");
  out_vx = fopen(FileName.c_str(),"w");
  //printf("passed second Filename overwirte\n");
  FILE*out_vy;
  FileName = "vy_t";
  FileName.append(bufferStep);
  FileName.append(".csv");
  out_vy = fopen(FileName.c_str(),"w");
  FILE*out_p;
  FileName = "p_t";
  FileName.append(bufferStep);
  FileName.append(".csv");
  out_p = fopen(FileName.c_str(),"w");
  FILE*out_Ma;
  FileName = "Ma_t";
  FileName.append(bufferStep);
  FileName.append(".csv");
  out_Ma = fopen(FileName.c_str(),"w");
  FILE*out_sos;
  FileName = "sos_t";
  FileName.append(bufferStep);
  FileName.append(".csv");
  out_sos = fopen(FileName.c_str(),"w");
  //printf("passed last Filename overwirte\n");
  for (int I = 0; I < Mx; I++)
    {
      for (int J = 0; J < My; J++)
	{
	  //printf("(I,J)=(%d,%d)\n",I,J); fflush(stdout);
	  fprintf(out_rho,"%8.7f\n", rho_local[I][J]);
	  fprintf(out_vx,"%8.7f\n", vx_local[I][J]);
	  fprintf(out_vy,"%8.7f\n", vy_local[I][J]);
	  fprintf(out_p,"%8.7f\n", p_local[I][J]);
	  fprintf(out_Ma,"%8.7f\n", sqrt(vx_local[I][J]*vx_local[I][J] + vy_local[I][J]*vy_local[I][J]) / sqrt(rsh*p_local[I][J]/rho_local[I][J]));
	  fprintf(out_sos,"%8.7f\n",sqrt(rsh*p_local[I][J]/rho_local[I][J]));
	}
    }
  fclose(out_rho);
  fclose(out_vx);
  fclose(out_vy);
  fclose(out_p);
  fclose(out_Ma);
  fclose(out_sos);
  printf("TRANSCRIPTION COMPLETE for step=%d\n\n", step);
}

int main()
{
  printf("\n\n***| ReorganizeData program has been launched |***\n\n");

  FILE*fsize = fopen("size_params.csv","r");
  int transfer_Mx;
  int transfer_My;
  fscanf(fsize, "%d,", &transfer_Mx);
  fscanf(fsize, "%d", &transfer_My);
  Mx = transfer_Mx;
  My = transfer_My;
  printf("---| Read in size_params.csv: Mx=%d, My=%d |---\n\n",Mx,My);
  
  
  for (int i = 0; i < 11; i++)
    {
      transcribe(i);
    }
  printf("\n\n***| Program execution successful, now exiting. |***\n\n");
  return 0;
}
