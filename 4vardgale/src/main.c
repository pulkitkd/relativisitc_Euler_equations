#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dg.h"
#include "dg1d.h"
#include "utils.h"


int main(int argc, char **argv) {

  // We want to print results only at specified times. The following code does
  // that, but it is very fragile. Bad inputs will break the code
  int step_output = 0;
  int skip_step = 5;
  if (argc > 1) {
    step_output = 1;
    skip_step = strtol(argv[1], &argv[1], 10);
  }

  UINT iter;
  REAL time;
  CELL *cell;
  FACE *face;

  GaussInit();
  cell = Init();
  face = InitFaces(cell);

  //   Result(cell); exit(0);
//printf("Debug: 1 \n");
  cfl = cfl / (2 * (PORD - 1) + 1);
  time = 0.0;
  iter = 0;
  // We set the error order to 4 orders more than PORD
  error_order = PORD + 4;

  printf("predictor_method = %d\n", predictor_method);
  printf("adaptation = %d\n", adaptation);

  printf("Beginning of iterations ...\n");
  if (step_output)
    printf("Writing output every %d iterations. \n", skip_step);

  // String for storing filename of the result
  while (time < finaltime) {
  if (step_output)
      if (iter % skip_step == 0)
        Result(cell, time, iter / skip_step, 0);
    // compute mesh velocity
    MeshVel(face);

    TimeStep(cell);
//printf("Debug: 2 \n");
    if (time + dt > finaltime)
      dt = finaltime - time;


    // Assemble residual
    Flux(cell, face, predictor_method);
//printf("Debug: 3 \n");
    // update solution
    Update(cell, face);
//printf("Debug: 4 \n");
    // Testing the merge now.
    // Dont adapt if this is last iteration
    if (adaptation == 1 && fabs(time+dt-finaltime)>1.0e-13)
      adapt(NCMAX, cell, NFMAX, face);
    // merge small cells, divide large

    // apply limiter
    ApplyLimiter(cell);
//printf("Debug: 5 \n");
    time += dt;
    ++iter;
    printf("%6d %14.6e %14.6e %14.6e %14.6e %10.3e %6d\n",
           iter,
           dt,
           time,
           dxmin,
           dxmax,
           cfl_f,
           NCMAX - free_cell_index
          );

    if (step_output)
      if (iter % skip_step == 0)
        Result(cell, time, iter/skip_step, 0);
    }

  Result(cell, time, iter, 1);
  int n_cell = 0;
  for (int i = 0; i < NCMAX; ++i)
    if (cell[i].active)
      ++n_cell;

  printf("Number of Active Cells = %d. \n", n_cell);
  printf("Number of Free Cells Remaining = %d. \n", NCMAX - free_cell_index);
  return 0;
}
