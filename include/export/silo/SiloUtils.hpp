#pragma once

// A collection of low-level utilities to help with silo file input/output.

/*-----------------------------------------------------------------------------
 * Purpose:     Impliment the create callback to initialize pmpio
 *              Will create the silo file and the 'first' directory (namespace)
 *              in it. The driver type (DB_PDB or DB_HDF5) is passed as user
 *              data; a void pointer to the driver determined in main.
 *-----------------------------------------------------------------------------
 */
void *CreateSiloFile(const char *fname, const char *nsname, void *userData)
{
    int driver = *((int*) userData);
    DBfile *siloFile = DBCreate(fname, DB_CLOBBER, DB_LOCAL, "pmpio testing", driver);
    if (siloFile)
    {
        DBMkDir(siloFile, nsname);
        DBSetDir(siloFile, nsname);
    }
    return (void *) siloFile;
}

/*-----------------------------------------------------------------------------
 * Purpose:     Impliment the open callback to initialize pmpio
 *              Will open the silo file and, for write, create the new
 *              directory or, for read, just cd into the right directory.
 *-----------------------------------------------------------------------------
 */
void *OpenSiloFile(const char *fname, const char *nsname, PMPIO_iomode_t ioMode,
    void *userData)
{
    DBfile *siloFile = DBOpen(fname, DB_UNKNOWN,
        ioMode == PMPIO_WRITE ? DB_APPEND : DB_READ);
    if (siloFile)
    {
        if (ioMode == PMPIO_WRITE)
            DBMkDir(siloFile, nsname);
        DBSetDir(siloFile, nsname);
    }
    return (void *) siloFile;
}

/*-----------------------------------------------------------------------------
 * Purpose:     Impliment the close callback for pmpio
 *-----------------------------------------------------------------------------
 */
void CloseSiloFile(void *file, void *userData)
{
    DBfile *siloFile = (DBfile *) file;
    if (siloFile)
        DBClose(siloFile);
}

// strdup isn't part of the C standard, so we can't rely on its existence.
// We keep our own handy.
char* strDup(const char* s)
{
  if (s == NULL)
    return NULL;
  char* dup = (char*)malloc(sizeof(char) * (strlen(s) + 1));
  strcpy(dup, s);
  return dup;
}

/*-----------------------------------------------------------------------------
 * Purpose:     Init silo file for pmpio
 *-----------------------------------------------------------------------------
 */
template<typename RealType>
void 
init_silo(diy::StaticAssigner& assigner,
	const std::string& directory, 
	int  timestep,
	RealType time,
	int numFiles,
	int mpiTag,
	DBfile *&silofile,
	PMPIO_baton_t *&baton,
	std::vector<std::vector<int> > &procsInGroup,
	std::vector<int> &nblockInRank)
{
  int driver = DB_HDF5;

  int nproc = 1, rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (numFiles == -1)
    numFiles = nproc;
  if (numFiles > nproc)
  {
    fprintf(stderr,"!!! ERROR: number of output Files > nproc\n");
    return;
  }

  baton = PMPIO_Init(numFiles, PMPIO_WRITE, MPI_COMM_WORLD, mpiTag, 
				    CreateSiloFile, 
				    OpenSiloFile, 
				    CloseSiloFile,
				    &driver);

  int groupRank = PMPIO_GroupRank(baton, rank);
  int rankInGroup = PMPIO_RankInGroup(baton, rank);

  procsInGroup.resize(numFiles);
  for (int i = 0; i < nproc; i++)
  {
    int groupRank = PMPIO_GroupRank(baton, i);
    procsInGroup[groupRank].push_back(i);
  }
  nblockInRank.resize(nproc);
  for (int i = 0; i < nproc; i++)
  {
    std::vector<int> gids;
    assigner.local_gids(i, gids);
    nblockInRank[i] = gids.size();
  }

  // Create the master directory if we need to.
  std::string masterDirName = directory;
  char dirname[1024];
  snprintf(dirname, 1024, "%s", masterDirName.c_str());
  masterDirName = dirname;
  
  if (rank == 0)
  {
    DIR* masterDir = opendir(masterDirName.c_str());
    if (masterDir == 0)
      mkdir((char*)masterDirName.c_str(), S_IRWXU | S_IRWXG);
    else
      closedir(masterDir);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  else
  {
    MPI_Barrier(MPI_COMM_WORLD);
  }

  char timestepdir[1024];
  snprintf(timestepdir, 1024, "%s/Timestep.%05d", masterDirName.c_str(), timestep);
  if (rankInGroup == 0)
  {
    DIR* groupDir = opendir(timestepdir);
    if (groupDir == 0)
      mkdir((char*)timestepdir, S_IRWXU | S_IRWXG);
    else
      closedir(groupDir);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  else
  {
    MPI_Barrier(MPI_COMM_WORLD);
  }

  char filename[2048];
  snprintf(filename, 2048, "%s/data_T%05d_G%05d.silo", timestepdir, timestep, groupRank);
  char domainame[1024];
  snprintf(domainame, 1024, "domain_%d", rankInGroup);
  silofile = (DBfile*)PMPIO_WaitForBaton(baton, filename, domainame);
  if(silofile == NULL)
  {
        fprintf(stderr, "Could not create Silo file!\n");
        return;
  }

}

/*-----------------------------------------------------------------------------
 * Purpose:     Write dumps.visit file for VISIT
 *-----------------------------------------------------------------------------
 */
void write_dumps(const std::string& directory, int  timestep)
{
  int nproc = 1, rank = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == nproc - 1)
  {	
      static bool s_summary_file_opened = false;
      std::string path = directory + "/dumps.visit";
      char temp[1024];
      snprintf(temp, 1024, "Timestep.%05d/summary.root", timestep);
      std::string file = temp;

      if (!s_summary_file_opened) {
         s_summary_file_opened = true;
         std::ofstream sfile(path.c_str(), std::ios::out);
         sfile <<  file << "\n";
         sfile.close(); 
      } else {
         std::ofstream sfile(path.c_str(), std::ios::app);
         sfile <<  file << "\n";
         sfile.close(); 
      }
  }
}



