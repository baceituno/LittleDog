function [plan, solvertime] = ral2017(quadruped, seed_plan, weights, goal_pos, use_symbolic)
% Runs a Mixed Integer Quadratic Program to plan footsteps on a Walking Robot.
% The structure of the desired footstep plan is indicated by the seed plan.
%
% @param seed_plan a FootstepPlan object which specifies the structure of the
%                  desired footstep plan. This seed_plan must contain a
%                  list of Footsteps and the region_order of matching length,
%                  but undetermined footstep positions and region assignments
%                  can be NaN. The first two footstep poses must not be NaN,
%                  since they correspond to the current positions of the feet
%
% @param weights a struct with fields 'goal', 'relative', and
%                'relative_final' describing the various contributions to
%                the cost function. These are described in detail in
%                hexapod.getFootstepOptimizationWeights()
%
% @param goal_pos a struct with fields 'Patai_Pie'.
%                 where goal_pos.Patai_Pie is the desired 6 DOF final pose
%                 of the ith foot
%
% @retval plan a FootstepPlan matching the number of footsteps, frame_id, etc.
%              of the seed plan, but with the footstep positions and region_order
%              replaced by the results of the MIQP
%
% @retval solvertime time taken by the solver (Gurobi by default) to generate the optimal footstep plan.

nsteps = length(seed_plan.footsteps);

t0 = tic();
planner = Quad_MixedIntegerFootstepPlanningProblem(quadruped, seed_plan);

planner.weights = weights;
%planner = planner.fixRotation();
planner = planner.addSinCosLinearEquality();
planner = planner.addReachabilityConstraints();
planner = planner.addGeometricConstraints();
%planner = planner.addTrimToFinalPoses();
%planner = planner.addQuadraticRelativeObjective();
planner = planner.addQuadraticGoalObjective(goal_pos, nsteps-3:nsteps, [1,1,1,1]);
planner = planner.addTerrainRegions([]);

fprintf(1, 'gurobi setup: %f\n', toc(t0));
[planner, solvertime, objval_nosymb] = planner.solve();
fprintf(1, 'gurobi total: %f\n', toc(t0));

plan = planner.getFootstepPlan();

end
