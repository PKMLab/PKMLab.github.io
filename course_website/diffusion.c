#include "build.h"
#include "ff.h"
#include "proto.h"
#include "diffusion.h"

main(int argc, char *argv[])
{
  int ans, i_mol, skip, i_atom, abs;

  int i_rel, j_rel, j_species;
  int i, j, k, i_atoms, i_species, i_phant, n_species, n_species_chk, 
      n_unit_cells, n_mols, n_atoms, n_mols1;
  int period_switch, vdw_switch, vdw_rad_switch, *excl_switch, init_switch,
      phant_switch, neigh_switch, opt_switch, press_switch, one_four_switch,
      nb_switch, con_switch, sim_switch, temp_switch, start_switch;

  /* Random generator */
  long i_ran, seed, idum;

  /* Array of geometrical molecular lengths. */
  double **mol_range, mol_max_length, xmax;

  /* Declare parameters used to build the initial state and box matrices*/
  int n_layers;
  double rho, aspect_xy, aspect_xz, layer_spacing;
  double **h, **h_inv, r_off, r_on, r2_on, r2_off, r_nl, r2_nl,
         volume;
  int kc;
  /* Declare and initialize force field parameters. */
  int n_mass = 0, n_lj = 0, n_e6 = 0, n_pot = 0, n_type = 0;
  mass_entry *mass = NULL;
  lj_entry *lj = NULL;
  e6_entry *e6 = NULL;
  pot_entry *pot = NULL;
  type_entry *atm_type = NULL;
  /* File names */
  char par_file[F_MAX], *mass_file, *lj_file, *e6_file, *pot_file, 
      **struct_file, *header_file, *config_file, *header_file_save, 
      *config_file_save;
  FILE *f_trajec, *f_com_pos, *f_output;
  /* Declare and initialize molecular structure parameters. */
  int *n_atoms_per_mol = NULL, *n_mols_per_species = NULL,
      *n_bonds_per_mol = NULL, **temp_atm_lab, **temp_atm_nbr,
      ***temp_atm_br, **temp_atm_type_i;
  double ***temp_atm_pos, **temp_atm_mass, **temp_atm_chg,
      **temp_atm_sigma, **temp_atm_eps, ***temp_atm_ord;
  char ***temp_atm_type;
  int **temp_bonds_1 = NULL, **temp_bonds_2 = NULL;
  double **temp_bonds_ord = NULL;
  int *n_mols_per_unit_cell, *mol_species;

  /* Declare some arrays used in set up force field */
  double *mol_mass, *mol_mass_true;
  int *mol_first_atm;

  /* Declare combination and exclusion arrays */
  int ***exclusions, **comb_pot;
  double ***comb_par;
  /* Declare arrays used in build */
  double **scaled_atom_coords, **mol_coords, **scaled_mol_coords,
      **atom_coords, **rel_atom_coords, **atom_vels, *atom_mass, 
      *sqrt_kT_over_m, *one_over_2m, **mol_vels;
  int *atom_rel, *atom_mol, *atom_type, *hot_atom;
  int hot_flag;

  /* Scan list arrays. */
  int **scan_atm_1, **scan_atm_2;

  /* Array of the chosen principal molecular axe. */
  int *principal;
  double **mol_inert_mt, ***mol_inert_axis;
  /* Monte Carlo parameters. */
  int i_cycle, n_cycles, n_trajec, n_block, n_blocks, n_inst;

/* Array for diffusion calculation */
   double **mol_coords_corr, **mol_coords_corr_0, **scaled_mol_coords_corr, delta;
   int n_msd, first_msd, i_data, water_diffusion;
   int pos_1, pos_2, record_size, i0, j0, j0max,
      i_msd, offset, *norm, n_fit, i_fit;
   double dr[NDIM], dr_p[NDIM], dum1, dum2, interval,
      *x, *y, *sig, a, b, chi2, q, siga, sigb, y_fit,
      diffus_x, sig_diffus_x,
      diffus_y, sig_diffus_y,
      diffus_z, sig_diffus_z,
      diffus_xy, sig_diffus_xy,
      diffus_xyz, sig_diffus_xyz,
      **r_displ2, **sig_r_displ2,
      *xy_displ2, *sig_xy_displ2,
      *xyz_displ2, *sig_xyz_displ2;

  double s_sep[NDIM], sep[NDIM], r2_sep;
    int j1, j2, rdfType, j_mol, n;
    char *type_i, *type_j;

  /* Get command line input. */
  if (argc != 2) {
    printf("Usage: %s params_file\n", argv[0]);
    exit(1);
  }
  strcpy(par_file, argv[1]);

  /* Read in input parameter file. */
  read_diffusion_params(par_file, &header_file, &config_file, 
      &mass_file, &pot_file, 
      &n_species_chk, &excl_switch, 
      &n_cycles, &n_block, &n_trajec, &seed, &delta, &n_msd, &first_msd, 
      &water_diffusion);

  printf("header_file = %s, config_file = %s\n",header_file,
          config_file);
  printf("mass_file = %s\n",mass_file);
  printf("pot_file = %s\n",pot_file);
  printf("n_species (from params file) = %d\n", n_species_chk);
  for (i_species = 0; i_species < n_species_chk; ++i_species)
       printf("i_specie = %d, excl_switch = %d\n", i_species,
          excl_switch[i_species]);

  printf("n_cycles = %d, \n", n_cycles);
  printf("n_trajec = %d\n", n_trajec);
  printf("seed = %ld\n",seed);
  printf("delta %g n_msd = %d first_msd %d\n", delta, n_msd, first_msd);

  /* Define additional variables. */
  if (n_cycles % n_block != 0) 
    error_exit("n_cycles must be a multiple of n_block!\n");
    n_blocks = n_cycles / n_block;

  /* Initialize random number generator. */
  i_ran = -seed;
  ran3(&i_ran);


  /* Read information about template molecules and atoms from header file */
  read_header_direct(header_file, &period_switch, &n_species,
      &n_atoms_per_mol,  &temp_atm_lab, &temp_atm_type, &temp_atm_nbr,
      &temp_atm_br, &temp_atm_ord,  &temp_atm_pos, &temp_atm_chg,
      &temp_atm_sigma, &temp_atm_eps,  &n_bonds_per_mol, &temp_bonds_1,
      &temp_bonds_2, &temp_bonds_ord,  &n_mols_per_species, &n_mols,
      &n_atoms, &mol_species);

  /* Check exclusion switch array. */
  if (n_species_chk != n_species) {
     printf("n_species_chk = %d, n_species = %d\n", n_species_chk, n_species);
     error_exit("Problem in the number of species.\n");
     }
 
  /* If the number of exclusion sites is greater than the number of atoms
     present in the template molecule, set the exclusion switch to -1, i.e.
     remove all intramolecular interactions. */ 
  for (i_species = 0; i_species < n_species; ++i_species)
    if (excl_switch[i_species] >= n_atoms_per_mol[i_species])
      excl_switch[i_species] = -1;
 
  /* Read FF masses file. */
  read_mass_params(mass_file, &mass, &n_mass);

  /* Set up force field data for masses. */
  set_up_force_field_pair(n_species, n_mols, n_atoms_per_mol, n_mass,
      mass, mol_species, n_bonds_per_mol, &temp_atm_mass, temp_atm_type,
      &mol_mass, &mol_mass_true, &mol_first_atm);

  /* Convert string tags to integer identifiers. */
  /* Create an integer identifier corresponding to template atom type. */ 
  convert_type(n_species, n_atoms_per_mol, temp_atm_type, &temp_atm_type_i, 
     &atm_type, &n_type);

  /* Set up scan list. */
  set_up_scan_list(n_species, n_atoms_per_mol, temp_atm_nbr, temp_atm_br,  
     &scan_atm_1, &scan_atm_2);


  /* Allocate memory for arrays used in the builder. */
  allocate_memory_diffusion(n_atoms, n_mols, n_species, 
      period_switch,  &scaled_atom_coords, &atom_rel, &atom_mol,
      &atom_mass, &atom_type, 
      &mol_coords, &scaled_mol_coords, 
      &rel_atom_coords, &h_inv, n_atoms_per_mol, &exclusions,
      &mol_inert_mt, &mol_inert_axis);


  /* Set up array of relative atom indices */
  relative_atoms_pair(n_mols, mol_species, mol_first_atm, n_atoms_per_mol,
      temp_atm_mass, temp_atm_type_i, atom_rel, atom_mol, atom_mass, atom_type);
      

  /* Set up the nonbonded exclusion array */
  exclusion_array(n_species, n_atoms_per_mol, excl_switch, temp_atm_nbr,  
    temp_atm_br, &exclusions);
 
   /* No of chromonic in the mixture */
  n_mols1 = n_mols_per_species[0];

   /* Allocate memory for and zero center of mass position and mean
      squared center of mass displacement arrays. */

       mol_coords_corr = allocate_2d_array(n_mols, 3, sizeof(double));
       mol_coords_corr_0 = allocate_2d_array(n_mols, 3, sizeof(double));
       r_displ2         = allocate_2d_array(n_msd+1, 3, sizeof(double));
       sig_r_displ2 = allocate_2d_array(n_msd+1, 3, sizeof(double));
       xy_displ2 = allocate_1d_array(n_msd+1, sizeof(double));
       sig_xy_displ2 = allocate_1d_array(n_msd+1, sizeof(double));
       xyz_displ2 = allocate_1d_array(n_msd+1, sizeof(double));
       sig_xyz_displ2 = allocate_1d_array(n_msd+1, sizeof(double));
       norm  = allocate_1d_array(n_msd+1, sizeof(int));
       scaled_mol_coords_corr = allocate_2d_array(n_mols, 3, sizeof(double));


  /* Open trajectory files. */
   f_trajec = gfopen(config_file, "r");
 /* Open output file to write unfolded center of mass */
   f_com_pos = gfopen("pro.com_pos_tmp", "w");

  /* Read header from the trajectory file. */
  fread(&n_atoms, sizeof(int), 1, f_trajec);
  fread(&n_cycles, sizeof(int), 1, f_trajec);
  fread(&n_trajec, sizeof(int), 1, f_trajec);

  pos_1 = ftell(f_trajec);

  printf("Header from %s\n", config_file);
  printf("n_atoms = %d\n", n_atoms);
  printf("n_cycles = %d\n", n_cycles);
  printf("n_trajec = %d\n", n_trajec);

  /* Compute number of instantaneous configuration in trajectory file. */
  if (n_cycles % n_trajec != 0) {
      printf("n_cycles must be a multiple of n_trajec\n");
      exit(1);
  }
  n_inst = n_cycles / n_trajec;

  /* Allocate memory for box and atom coordinates. */
  h = allocate_2d_array(3, 3, sizeof(double));
  atom_coords = allocate_2d_array(8*n_atoms, 3, sizeof(double));

/* calculate unfolded molecular center of mass positions */

  /* Position trajectory file pointer at beginning of first data set. */
      fseek(f_trajec, pos_1, 0);


   /* Read in first data set. */
      read_positions_direct(f_trajec, h, n_atoms, atom_coords);

   /* If params.period_switch == 1, calculate quantities that depend
      on box dimensions, calculate scaled atomic coordinates, and apply
      periodic boundary conditions. */
   if (period_switch) {
      box_dimensions(h, h_inv, period_switch, 0.0, 0.0);
      scaled_atomic_coords(n_atoms, h_inv, atom_coords, scaled_atom_coords);
      periodic_boundary_conditions(n_atoms, h, scaled_atom_coords,
                                      atom_coords);
   }

   /* Calculate center of mass positions. */
      center_of_mass_positions(period_switch, n_mols, n_species,
      n_atoms_per_mol,  mol_species, atom_mass, mol_first_atm,
      atom_coords, scaled_atom_coords,  scan_atm_1, scan_atm_2,
      mol_mass, h, h_inv, mol_coords, scaled_mol_coords,  rel_atom_coords);

  /* If period_switch == 1, initialize corrected center of mass
      positions. */
   if (period_switch)
      for (i_mol = 0; i_mol < n_mols; ++i_mol) {
          for ( k = 0; k < NDIM; ++k) {
         scaled_mol_coords_corr[i_mol][k] = scaled_mol_coords[i_mol][k];
         mol_coords_corr[i_mol][k] = mol_coords[i_mol][k];
        }
      }

   /* If period_switch == 1, write out corrected center of mass
      positions. Otherwise, write out "raw" center of mass positions. */
  for (i_mol = 0; i_mol  < n_mols; ++i_mol) {
   if (period_switch)
      fwrite(mol_coords_corr[i_mol], sizeof(double), 3, f_com_pos);
   else
      fwrite(mol_coords[i_mol], sizeof(double), 3, f_com_pos);
 }
   /* Loop over data sets. */
   for (i_data = 1; i_data < n_inst; ++i_data) {

      /* Read in atomic positions. */
      read_positions_direct(f_trajec, h, n_atoms, atom_coords);

      /* If params.period_switch == 1, calculate quantities that depend
         on box dimensions, calculate scaled atomic coordinates, and apply
         periodic boundary conditions. */
      if (period_switch) {
      box_dimensions(h, h_inv, period_switch, 0.0, 0.0);
      scaled_atomic_coords(n_atoms, h_inv, atom_coords, scaled_atom_coords);
      periodic_boundary_conditions(n_atoms, h, scaled_atom_coords,
                                      atom_coords);
      }

      /* Calculate center of mass positions. */
      center_of_mass_positions(period_switch, n_mols, n_species,
      n_atoms_per_mol,  mol_species, atom_mass, mol_first_atm,
      atom_coords, scaled_atom_coords,  scan_atm_1, scan_atm_2,
      mol_mass, h, h_inv, mol_coords, scaled_mol_coords,  rel_atom_coords);
 
      /* If period_switch == 1, calculate corrected scaled
         center of mass positions. */

      if (period_switch)
         for (i_mol = 0; i_mol < n_mols; ++i_mol) {
             for (k = 0; k <NDIM; ++k){
            sep[k] = scaled_mol_coords[i_mol][k] - scaled_mol_coords_corr[i_mol][k];
            sep[k] -= NINT(sep[k]);
            scaled_mol_coords_corr[i_mol][k] += sep[k];
            mol_coords_corr[i_mol][k] =
               h[k][k] * scaled_mol_coords_corr[i_mol][k];
           }
         }

      /* If period_switch == 1, write out corrected center of mass
         positions. Otherwise, write out "raw" center of mass positions. */
    for ( i_mol = 0; i_mol < n_mols; ++i_mol) {
      if (period_switch)
         fwrite(mol_coords_corr[i_mol], sizeof(double), 3, f_com_pos);
      else
         fwrite(mol_coords[i_mol], sizeof(double), 3, f_com_pos);
    }
   }

   /* Close input and output files. */
   fclose(f_trajec);
   fclose(f_com_pos);

  /* Calculate mean Squared displacement from the unfolded coordinate */

  /* Determine size (in bytes) of records. reopen the unfolded center of mass file */

   f_com_pos = gfopen("pro.com_pos_tmp", "r");

   pos_1 = ftell(f_com_pos);
   for (i_mol = 0; i_mol < n_mols; ++i_mol)
   fread(mol_coords_corr[i_mol], sizeof(double), 3, f_com_pos);
   pos_2 = ftell(f_com_pos);
   record_size = pos_2 - pos_1;
   fclose(f_com_pos);

   /* Exit if n_msd > inst - 1. */
   if (n_msd > n_inst - 1)
      error_exit("n_msd > n_inst - 1 in mean_squared_com_displacement.");

   /* Reopen input file. */
    f_com_pos = gfopen("pro.com_pos_tmp", "r");

   /* Add contributions to mean squared center of mass displacement. */
   for (i0 = 0; i0 < n_inst; ++i0) {

      /* Calculate file position offset for initial data set, and read data
         set starting at that point in the data file. */
      offset = i0 * record_size;
      fseek(f_com_pos, offset, 0);
      for (i_mol = 0; i_mol < n_mols; ++i_mol)
         fread(mol_coords_corr_0[i_mol], sizeof(double), 3, f_com_pos);

      /* Loop over time offsets. */
      j0max = MIN(n_inst - 1, i0 + n_msd);
      for (j0 = i0; j0 <= j0max; ++j0) {

         /* Calculate file position offset for final data set, and read data
            set starting at that point in the data file. */
         offset = j0 * record_size;
         fseek(f_com_pos, offset, 0);
         for (i_mol = 0; i_mol < n_mols; ++i_mol)
           fread(mol_coords_corr[i_mol], sizeof(double), 3, f_com_pos);

/*    Now chek for the species */
      if(water_diffusion == 0) n_mols = n_mols1;
        else n_mols = n_mols - n_mols1;

      /* Add contributions to mean squared center of mass displacement. */
         i_msd = j0 - i0;
         for (i_mol = 0; i_mol < n_mols; ++i_mol) {
             for (k = 0; k < NDIM; ++k){
                dr[k] = mol_coords_corr[i_mol][k] - mol_coords_corr_0[i_mol][k];
                dr_p[k] = dr[k] * h[k][k];
                r_displ2[i_msd][k] += SQR(dr_p[k]);
             }
         }
         ++norm[i_msd];
      }
   }

   /* Close input file. */
   fclose(f_com_pos);

  /* Normalize mean squared center of mass displacement, calculate
      uncertainties, and write results to output file. Some remarks are
      in order regarding the calculation of uncertainties. The uncertainty
      in <[x(t)-x(0)]^2> is assumed to be equal to the standard deviation
      in [x(t)-x(0)]^2 (which is equal to sqrt(2) * <[x(t)-x(0)]^2> for a
      Gaussian distribution) divided by the square root of the number of
      independent measurements of <[x(t)-x(0)]^2>. Following Allen and
      Tildesley, the number of independent measurements is taken to be
      the total sampling interval divided by twice the correlation
      time for [x(t)-x(0)]^2 , which is estimated to be equal to t (valid
      for large t). We have further assumed that each molecule makes an
      independent contribution to the mean squared displacement. */

   interval = delta * n_trajec;

   f_output = gfopen("pro.msd", "w");

   for (i_msd = 0; i_msd <= n_msd; ++i_msd) {
      dum1 = 1.0 / (n_mols * norm[i_msd]);

      for(k = 0; k < NDIM; ++k)
          r_displ2[i_msd][k] *= dum1;

      xy_displ2[i_msd] = r_displ2[i_msd][0] + r_displ2[i_msd][1];
      xyz_displ2[i_msd] = r_displ2[i_msd][0] + r_displ2[i_msd][1] + r_displ2[i_msd][2];
      dum2 = sqrt((4.0 * i_msd) / (n_mols * (n_inst - i_msd)));

      for(k = 0; k < NDIM; ++k)
          sig_r_displ2[i_msd][k] = dum2 * r_displ2[i_msd][k];

      sig_xy_displ2[i_msd] = sqrt(SQR(sig_r_displ2[i_msd][0])
         + SQR(sig_r_displ2[i_msd][1]));

      sig_xyz_displ2[i_msd] = sqrt(SQR(sig_r_displ2[i_msd][0])
         + SQR(sig_r_displ2[i_msd][1]) + SQR(sig_r_displ2[i_msd][2]));

    fprintf(f_output, "%g %g %g %g %g %g %g %g %g %g %g\n",
         i_msd * interval,
         r_displ2[i_msd][0], sig_r_displ2[i_msd][0],
         r_displ2[i_msd][1], sig_r_displ2[i_msd][1],
         r_displ2[i_msd][2], sig_r_displ2[i_msd][2],
         xy_displ2[i_msd], sig_xy_displ2[i_msd],
         xyz_displ2[i_msd], sig_xyz_displ2[i_msd]);
   }
   fclose(f_output);

   /* Calculate diffusion constants from linear fits to the mean square
      displacements. */
   x = dvector(1, n_msd + 1);
   y = dvector(1, n_msd + 1);
   sig = dvector(1, n_msd + 1);
   n_fit = n_msd - first_msd + 1;
   for (i_fit = 1; i_fit <= n_fit; ++i_fit) {
      i_msd = i_fit + first_msd - 1;
      x[i_fit] = i_msd * interval;
      y[i_fit] = r_displ2[i_msd][0];
      sig[i_fit] = sig_r_displ2[i_msd][0];
   }
   fit(x, y, n_fit, sig, 1, &a, &b, &siga, &sigb, &chi2, &q);
   diffus_x = b / 2.0;
   sig_diffus_x = sigb / 2.0;
   f_output = fopen("pro.msd_fit_x", "w");
   for (i_fit = 1; i_fit <= n_fit; ++i_fit) {
      y_fit = a + b * x[i_fit];
      fprintf(f_output, "%g %g %g %g\n", x[i_fit], y[i_fit], sig[i_fit], y_fit);
   }
   fclose(f_output);
   for (i_fit = 1; i_fit <= n_fit; ++i_fit) {
      i_msd = i_fit + first_msd - 1;
      x[i_fit] = i_msd * interval;
      y[i_fit] = r_displ2[i_msd][1];
      sig[i_fit] = sig_r_displ2[i_msd][1];
   }
   fit(x, y, n_fit, sig, 1, &a, &b, &siga, &sigb, &chi2, &q);
   diffus_y = b / 2.0;
   sig_diffus_y = sigb / 2.0;
   f_output = fopen("pro.msd_fit_y", "w");
   for (i_fit = 1; i_fit <= n_fit; ++i_fit) {
      y_fit = a + b * x[i_fit];
      fprintf(f_output, "%g %g %g %g\n", x[i_fit], y[i_fit], sig[i_fit], y_fit);
   }
   fclose(f_output);
   for (i_fit = 1; i_fit <= n_fit; ++i_fit) {
      i_msd = i_fit + first_msd - 1;
      x[i_fit] = i_msd * interval;
      y[i_fit] = r_displ2[i_msd][2];
      sig[i_fit] = sig_r_displ2[i_msd][2];
   }
   fit(x, y, n_fit, sig, 1, &a, &b, &siga, &sigb, &chi2, &q);
   diffus_z = b / 2.0;
   sig_diffus_z = sigb / 2.0;
   f_output = fopen("pro.msd_fit_z", "w");
   for (i_fit = 1; i_fit <= n_fit; ++i_fit) {
      y_fit = a + b * x[i_fit];
      fprintf(f_output, "%g %g %g %g\n", x[i_fit], y[i_fit], sig[i_fit], y_fit);
   }
   fclose(f_output);
   for (i_fit = 1; i_fit <= n_fit; ++i_fit) {
      i_msd = i_fit + first_msd - 1;
      x[i_fit] = i_msd * interval;
      y[i_fit] = xy_displ2[i_msd];
      sig[i_fit] = sig_xy_displ2[i_msd];
   }
   fit(x, y, n_fit, sig, 1, &a, &b, &siga, &sigb, &chi2, &q);
   diffus_xy = b / 4.0;
   sig_diffus_xy = sigb / 4.0;
   f_output = fopen("pro.msd_fit_xy", "w");
   for (i_fit = 1; i_fit <= n_fit; ++i_fit) {
      y_fit = a + b * x[i_fit];
      fprintf(f_output, "%g %g %g %g\n", x[i_fit], y[i_fit], sig[i_fit], y_fit);
   }
   fclose(f_output);
   for (i_fit = 1; i_fit <= n_fit; ++i_fit) {
      i_msd = i_fit + first_msd - 1;
      x[i_fit] = i_msd * interval;
      y[i_fit] = xyz_displ2[i_msd];
      sig[i_fit] = sig_xyz_displ2[i_msd];
   }
   fit(x, y, n_fit, sig, 1, &a, &b, &siga, &sigb, &chi2, &q);
   diffus_xyz = b / 6.0;
   sig_diffus_xyz = sigb / 6.0;
   f_output = fopen("pro.msd_fit_xyz", "w");
   for (i_fit = 1; i_fit <= n_fit; ++i_fit) {
      y_fit = a + b * x[i_fit];
      fprintf(f_output, "%g %g %g %g\n", x[i_fit], y[i_fit], sig[i_fit], y_fit);
   }
   fclose(f_output);

   /* Write diffusion constants to output file. */
   dum1 = 1.0e-4;
   dum2 = delta;
   f_output = fopen("pro.diffusion", "w");
   fprintf(f_output, "diffus_x = %g +/- %g angstrom^2 / ps\n",
      diffus_x, sig_diffus_x);
   fprintf(f_output, "         = %g +/- %g cm^2 / s\n",
      dum1 * diffus_x, dum1 * sig_diffus_x);
   fprintf(f_output, "         = %g +/- %g angstrom^2 / MD step\n",
      dum2 * diffus_x, dum2 * sig_diffus_x);
   fprintf(f_output, "diffus_y = %g +/- %g angstrom^2 / ps\n",
      diffus_y, sig_diffus_y);
   fprintf(f_output, "         = %g +/- %g cm^2 / s\n",
      dum1 * diffus_y, dum1 * sig_diffus_y);
   fprintf(f_output, "         = %g +/- %g angstrom^2 / MD step\n",
      dum2 * diffus_y, dum2 * sig_diffus_y);
   fprintf(f_output, "diffus_z = %g +/- %g angstrom^2 / ps\n",
      diffus_z, sig_diffus_z);
   fprintf(f_output, "         = %g +/- %g cm^2 / s\n",
      dum1 * diffus_z, dum1 * sig_diffus_z);
   fprintf(f_output, "         = %g +/- %g angstrom^2 / MD step\n",
      dum2 * diffus_z, dum2 * sig_diffus_z);
   fprintf(f_output, "diffus_xy = %g +/- %g angstrom^2 / ps\n",
      diffus_xy, sig_diffus_xy);
   fprintf(f_output, "          = %g +/- %g cm^2 / s\n",
      dum1 * diffus_xy, dum1 * sig_diffus_xy);
   fprintf(f_output, "          = %g +/- %g angstrom^2 / MD step\n",
        dum2 * diffus_xy, dum2 * sig_diffus_xy);
   fprintf(f_output, "diffus_xyz = %g +/- %g angstrom^2 / ps\n",
      diffus_xyz, sig_diffus_xyz);
   fprintf(f_output, "           = %g +/- %g cm^2 / s\n",
      dum1 * diffus_xyz, dum1 * sig_diffus_xyz);
   fprintf(f_output, "           = %g +/- %g angstrom^2 / MD step\n",
      dum2 * diffus_xyz, dum2 * sig_diffus_xyz);
   fclose(f_output);

   /* Free up memory. */
   free(mol_coords_corr);
   free(mol_coords_corr_0);
   free(r_displ2);
   free(sig_r_displ2);
   free(xy_displ2);
   free(sig_xy_displ2);
   free(xyz_displ2);
   free(sig_xyz_displ2);
   free(norm);
   free_dvector(x, 1, n_fit);
   free_dvector(y, 1, n_fit);
   free_dvector(sig, 1, n_fit);

  exit(0);

}  /* end of main. */

/****************************************************************************/
/* Allocate memory for static arrays */
void allocate_memory_diffusion(int n_atoms, int n_mols, int n_species,
int period_switch, double ***scaled_atom_coords, int **atom_rel,
int **atom_mol, double **atom_mass, int **atom_type, 
double ***mol_coords, 
double ***scaled_mol_coords, double ***rel_atom_coords, 
double ***h_inv, int *n_atoms_per_mol,int ****exclusions, 
double ***mol_inert_mt, double ****mol_inert_axis)

{
  int i, i_rel, i_species, n_atoms_max;

  /* Allocate memory for arrays of atomic properties. */
  if (period_switch)
    *scaled_atom_coords = allocate_2d_array(n_atoms,
        3, sizeof(double));
  *rel_atom_coords = allocate_2d_array(n_atoms, 3, sizeof(double));

  /* Allocate memory for arrays of atomic properties. */
  *atom_rel = allocate_1d_array(n_atoms, sizeof(int));
  *atom_mol = allocate_1d_array(n_atoms, sizeof(int));
  *atom_mass = allocate_1d_array(n_atoms, sizeof(double));
  *atom_type = allocate_1d_array(n_atoms, sizeof(int));

  /* Allocate memory for arrays of molecular properties. */
  *mol_coords = allocate_2d_array(n_mols, 3, sizeof(double));
  if (period_switch)
    *scaled_mol_coords = allocate_2d_array(n_mols, 3,
        sizeof(double));

  /* Allocate memory for the box matrices */
  /* For now, orthorombic box */
  *h_inv = allocate_2d_array(3, 3, sizeof(double));

  /* Allocate memory for exclusion array. */
  *exclusions = allocate_1d_array(n_species, sizeof(int**));
  for (i_species = 0; i_species < n_species; ++i_species) {
    (*exclusions)[i_species] = allocate_1d_array(n_atoms_per_mol[i_species],
        sizeof(int*));
    for (i_rel = 0; i_rel < n_atoms_per_mol[i_species];
        ++i_rel) {
      (*exclusions)[i_species][i_rel] = allocate_1d_array(n_atoms_per_mol[i_species],
          sizeof(int));
    }
  }
   *mol_inert_mt = allocate_2d_array(n_mols, 3, sizeof(double));
   *mol_inert_axis = allocate_3d_array(n_mols, 3, 3, sizeof(double));
}

/****************************************************************************/
/* Set up array of relative atom indices. */
void relative_atoms_pair(int n_mols, int *mol_species, int *mol_first_atm,
   int *n_atoms_per_mol, double **temp_atm_mass, int **temp_atm_type_i,
   int *atom_rel, int *atom_mol, double *atom_mass, int *atom_type)
{
  int i_species, skip, i, i_mol, i_rel;

  /* Set up array of atomic properties. */
  for (i_mol = 0; i_mol < n_mols; ++i_mol) {
    i_species = mol_species[i_mol];
    skip = mol_first_atm[i_mol];
    for (i_rel = 0; i_rel < n_atoms_per_mol[i_species]; ++i_rel) {
      i = skip + i_rel;
      atom_rel[i] = i_rel;
      atom_mol[i] = i_mol;
      atom_mass[i] = temp_atm_mass[i_species][i_rel];
      atom_type[i] = temp_atm_type_i[i_species][i_rel];
    }
  }

}
