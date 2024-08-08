var documenterSearchIndex = {"docs":
[{"location":"reference/api/#API-Reference","page":"API Reference","title":"API Reference","text":"","category":"section"},{"location":"reference/api/","page":"API Reference","title":"API Reference","text":"Modules = [Jessamine]\nOrder   = [:module, :constant, :type, :function, :macro,]","category":"page"},{"location":"reference/api/#Jessamine.AbstractGeneOp","page":"API Reference","title":"Jessamine.AbstractGeneOp","text":"An operation to be performed as part of the function of a gene.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.AbstractGenome","page":"API Reference","title":"Jessamine.AbstractGenome","text":"Abstract base type for genomes.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.AbstractLinearModelResult","page":"API Reference","title":"Jessamine.AbstractLinearModelResult","text":"Abstract base type for the results of fitting linear models.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.AbstractMachineSpec","page":"API Reference","title":"Jessamine.AbstractMachineSpec","text":"Abstract base type for machine specs\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.AbstractModelResult","page":"API Reference","title":"Jessamine.AbstractModelResult","text":"Abstract base type for the results of fitting a model to data.  This object can then be used for predictions.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.AbstractPopulationCondition","page":"API Reference","title":"Jessamine.AbstractPopulationCondition","text":"Base type of conditions for a population.  These indicate whether the population is a work in progress, or why its evolution was stopped.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.AbstractSolverSpec","page":"API Reference","title":"Jessamine.AbstractSolverSpec","text":"Abstract base type for solver specs, used for solving for genome parameter values.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.Add","page":"API Reference","title":"Jessamine.Add","text":"Add operands.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.Agent","page":"API Reference","title":"Jessamine.Agent","text":"An agent has a rating, a genome, a parameter vector, and an extra bit of data that depends on how rating is done.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.BasicLinearModelResult","page":"API Reference","title":"Jessamine.BasicLinearModelResult","text":"A BasicLinearModelResult has a vector of coefficients plust a bias or intercept.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.DefaultSolverSpec","page":"API Reference","title":"Jessamine.DefaultSolverSpec","text":"Use this if you want all default solving procedures.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.EvolutionSpec","page":"API Reference","title":"Jessamine.EvolutionSpec","text":"Parameters needed to run random_initial_population and evolution_loop.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.FzAnd","page":"API Reference","title":"Jessamine.FzAnd","text":"Return fuzzy AND of the operands\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.FzNand","page":"API Reference","title":"Jessamine.FzNand","text":"Return fuzzy NAND of the operands.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.FzNor","page":"API Reference","title":"Jessamine.FzNor","text":"Return fuzzy NOR of the operands.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.FzOr","page":"API Reference","title":"Jessamine.FzOr","text":"Return fuzzy OR of the operands.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.Genome","page":"API Reference","title":"Jessamine.Genome","text":"A vector of blocks of instructions.  During a time step, each instruction is evaluated on the current work space.  For each j, the results of the instructions in block j are collected and added, and this sum is used as the value of slot j in the next work space.  \n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.GenomeSpec","page":"API Reference","title":"Jessamine.GenomeSpec","text":"A collection of parameters specifying the genome architecture and mutation processes.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.LinearModelMachineSpec","page":"API Reference","title":"Jessamine.LinearModelMachineSpec","text":"Wrap an MLJ Model that uses a linear combination of input columns to produce its predictions.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.LinearModelMachineSpec-Tuple{MLJModelInterface.Model, Any, Any}","page":"API Reference","title":"Jessamine.LinearModelMachineSpec","text":"LinearModelMachineSpec(model::Model, lambda_parameter, lambda_operand)\n\nConstruct a LinearModelMachineSpec using model.lambda for lambda_model.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.Maximum","page":"API Reference","title":"Jessamine.Maximum","text":"Return the maximum of the operands.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.Minimum","page":"API Reference","title":"Jessamine.Minimum","text":"Return the minimum of the operands.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.Multiply","page":"API Reference","title":"Jessamine.Multiply","text":"Multiply operands.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.MutationDist","page":"API Reference","title":"Jessamine.MutationDist","text":"A collection of distribution objects built from a MutationSpec and used to produce random numbers determining mutations.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.MutationDist-Tuple{MutationSpec, Int64}","page":"API Reference","title":"Jessamine.MutationDist","text":"MutationDist(m_spec::MutationSpec, src_index_max::Int)\n\nUse m_spec to build a MutationDist. Specify that an operand index must be <= src_index_max.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.MutationSpec","page":"API Reference","title":"Jessamine.MutationSpec","text":"A collection of parameters specifying the probabilities of various mutations.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.Population","page":"API Reference","title":"Jessamine.Population","text":"A population is a single generation of living Agents.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.RDDState","page":"API Reference","title":"Jessamine.RDDState","text":"State structure corresponding to ongoing iteration of RandomDuplicateDelete\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.RandomDuplicateDelete","page":"API Reference","title":"Jessamine.RandomDuplicateDelete","text":"Wrap an iterator source to produce an iterator that will randomly include some elements yielded by the source twice in a row, and some not at all.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.ReciprocalAdd","page":"API Reference","title":"Jessamine.ReciprocalAdd","text":"Add operands and return the reciprocal.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.ReciprocalMultiply","page":"API Reference","title":"Jessamine.ReciprocalMultiply","text":"Multiply operands and return the reciprocal.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.ReciprocalSubtract","page":"API Reference","title":"Jessamine.ReciprocalSubtract","text":"Subtract operands and return the reciprocal.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.SelectionDist","page":"API Reference","title":"Jessamine.SelectionDist","text":"Distribution object for tournament selection.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.SelectionSpec","page":"API Reference","title":"Jessamine.SelectionSpec","text":"Parameters for tournament selection.\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.SignAdd","page":"API Reference","title":"Jessamine.SignAdd","text":"Return the sign of the sum of the operands\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.Subtract","page":"API Reference","title":"Jessamine.Subtract","text":"Subtract LISP-style: 0 + x[1] - x[2] - x[3]...\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Jessamine.UnaryComposition","page":"API Reference","title":"Jessamine.UnaryComposition","text":"Do a multi-arity operation and apply a unary operation\n\n\n\n\n\n","category":"type"},{"location":"reference/api/#Base.isless-Tuple{Agent, Agent}","page":"API Reference","title":"Base.isless","text":"isless(x::Agent, y::Agent)\n\nCompare the rating field of the two agents. This definition allows a population to be sorted.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.area_above_curve-Tuple{Any, Any}","page":"API Reference","title":"Jessamine.area_above_curve","text":"area_above_curve(y_hat, y_ref)\n\nReturn the area above the ROC curve for the predictions y_hat and actual reference values y_ref. This is 1 minus the area under the curve. The vector y_hat should consist of Bernoulli distributions, as returned by a LogisticClassifier, for example. The vector y should consist of values from those distributions.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.describe_condition-Tuple{AbstractPopulationCondition}","page":"API Reference","title":"Jessamine.describe_condition","text":"describe_condition(::AbstractPopulationCondition)\n\nReturn a short string describing the population condition.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.do_simplification-Tuple{AbstractPopulationCondition}","page":"API Reference","title":"Jessamine.do_simplification","text":"do_simplification(::AbstractPopulationCondition)\n\nReturn whether to perform the simplification epoch. Utility function.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.eval_time_step_symbolic-Tuple{GenomeSpec, Genome}","page":"API Reference","title":"Jessamine.eval_time_step_symbolic","text":"eval_time_step_symbolic(g_spec, genome; output_sym, scratch_sym, parameter_sym, input_sym)\n\nReturn symbolic objects representing a single time step of genome. The result is a named tuple with fields z, t, p, x for the Symbolics objects for those variables; c and c_next for the current and future cell state in symbolic form.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.extend_if_singleton-Tuple{AbstractVector, Int64}","page":"API Reference","title":"Jessamine.extend_if_singleton","text":"extend_if_singleton(v::AbstractVector, m::Int)\n\nIf v is a singleton, as in v = [x], return [x, x, ...], that is, a vector of m copies of x. Otherwise, assert that length(v) == m and return v. So the result is always a vector of length m.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.generation_size-Tuple{SelectionSpec}","page":"API Reference","title":"Jessamine.generation_size","text":"generation_size(s_spec::SelectionSpec)\n\nReturn the size of a generation.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.genome_complexity-Tuple{AbstractMachineSpec, GenomeSpec, Genome}","page":"API Reference","title":"Jessamine.genome_complexity","text":"genome_complexity(mn_spec, g_spec, genome)\n\nNumerical complexity of genome.  The default implementation returns mn_spec.lambda_operand times the number of operands in the genome.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.genome_parameter_complexity-Tuple{AbstractMachineSpec, AbstractVector}","page":"API Reference","title":"Jessamine.genome_parameter_complexity","text":"genome_parameter_complexity(mn_spec, p)\n\nNumerical complexity of a parameter vector.  The default implementation returns mn_spec.lambda_parameter times the squared L2 norm of p.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.least_squares_ridge-Tuple{AbstractVector{<:AbstractVector}, AbstractVector, Number, GenomeSpec, AbstractGenome, AbstractVector}","page":"API Reference","title":"Jessamine.least_squares_ridge","text":"least_squares_ridge(xs, y, lambda, g_spec, genome, parameter)\n\nCompute ridge regression. Assume xs is an array of columns as predictors and y is a column of target values. Apply run_genome using parameter as the parameter vector and xs as the inputs. Gather a column of 1s and the output columns as a matrix X. The prediction values are y_hat = X * b, where b is a column of (unknown) coefficients. Solve for the b that minimizes norm(y - y_hat)^2 + lambda * norm(b)^2. If all goes well, return (norm(y - y_hat), b). Otherwise return nothing.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.least_squares_ridge_grow_and_rate","page":"API Reference","title":"Jessamine.least_squares_ridge_grow_and_rate","text":"least_squares_ridge_grow_and_rate(xs, y, lambda_b, lambda_p, lambda_op, g_spec, genome, p_init = zeros(...))\n\nSolve for the parameter vector p that minimzes norm(y - y_hat)^2 + lambda_b * norm(b)^2 + lambda_p * norm(p)^2 + lambda_op R, where y_hat and b are found using least_squares_ridge. R is the total number of operands across all instructions in genome.\n\nThe solver starts with p_init for the initial value of p. If p_init is nothing or not given, a vector of zeros is used.\n\nIf all goes well, return an Agent, whose genome is genome, whose parameter is the best p, and whose extra is a BasicLinearModelResult with coefficient vector b.\n\nOtherwise, return nothing.\n\n\n\n\n\n","category":"function"},{"location":"reference/api/#Jessamine.machine_complexity-Tuple{AbstractMachineSpec, Any}","page":"API Reference","title":"Jessamine.machine_complexity","text":"machine_complexity(mn_spec, m; kw_args...)\n\nNumerical complexity of a machine.\n\nThe default is to return 0.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.machine_complexity-Tuple{LinearModelMachineSpec, Any}","page":"API Reference","title":"Jessamine.machine_complexity","text":"machine_complexity(mn_spec, m)\n\nNumerical complexity of a machine.\n\nReturn the sum of squares of MLJ.fitted_params(m) multiplied by mn_spec.lambda_model.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.machine_fit!-Tuple{AbstractMachineSpec, Any}","page":"API Reference","title":"Jessamine.machine_fit!","text":"machine_fit!(mn_spec, m; kw_args...)\n\nThe default implementation returns MLJ.fit!(m; kw_args...).\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.machine_grow_and_rate","page":"API Reference","title":"Jessamine.machine_grow_and_rate","text":"machine_grow_and_rate(xs, y, g_spec, genome, mn_spec, sol_spec, w = nothing)\n\nSolve for the parameter values that produce the best machine for predicting y from the outputs of the genome applied to xs. A non-nothing value of w is passed to the machine as a weight vector. The resulting parameter vector is stored as the parameter field in the returned Agent.  The resulting machine is wrapped in a MachineResult and stored in the extra field of the Agent.\n\nIf any exception is thrown during the solving process, the exception is suppressed, and nothing is returned.\n\n\n\n\n\n","category":"function"},{"location":"reference/api/#Jessamine.machine_init-Tuple{AbstractMachineSpec, Any, Any}","page":"API Reference","title":"Jessamine.machine_init","text":"machine_init(mn_spec, X, y; kw_args...)\n\nConstruct a machine with predictor table X and target y.\n\nThe default implementation calls MLJ.machine(mn_spec.machine, X, y; kw_args...)\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.machine_predict-Tuple{AbstractMachineSpec, Any, Any}","page":"API Reference","title":"Jessamine.machine_predict","text":"machine_predict(mn_spec, m, X; kw_args...)\n\nProduce a prediction y_hat. The default implementation returns MLJ.predict(m, X; kw_args...).\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.machine_spec_type-Tuple{Type{MLJLinearModels.LassoRegressor}}","page":"API Reference","title":"Jessamine.machine_spec_type","text":"machine_spec_type(::Type{LassoRegressor})\n\nReturn LinearModelMachineSpec\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.machine_spec_type-Tuple{Type{MLJLinearModels.LogisticClassifier}}","page":"API Reference","title":"Jessamine.machine_spec_type","text":"machine_spec_type(::Type{LogisticClassifier})\n\nReturn LinearModelMachineSpec\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.machine_spec_type-Tuple{Type{MLJLinearModels.RidgeRegressor}}","page":"API Reference","title":"Jessamine.machine_spec_type","text":"machine_spec_type(::Type{RidgeRegressor})\n\nReturn LinearModelMachineSpec\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.machine_spec_type-Tuple{Type}","page":"API Reference","title":"Jessamine.machine_spec_type","text":"machine_spec_type(t_MLJ_model::Type)::Type{<:AbstractMachineSpec}\n\nReturn the Jessamine type (subtype of AbstractMachineSpec) corresponding to the given MLJ model type. The default is to return BasicModelMachineSpec\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.model_predict","page":"API Reference","title":"Jessamine.model_predict","text":"model_predict(mr::AbstractModelResult, X; kw_args...)\n\nGiven the table X where the columns are inputs (predictors) and the rows are points, use mr to predict the target value for all the points.\n\n\n\n\n\n","category":"function"},{"location":"reference/api/#Jessamine.model_predict-Tuple{AbstractLinearModelResult, AbstractVector{<:AbstractVector}}","page":"API Reference","title":"Jessamine.model_predict","text":"model_predict(mr::AbstractLinearModelResult, xs::AbstractVector{<:AbstractVector}; kw_args...)\n\nStack the columns xs into a matrix X and call model_predict. The kw_args are splatted in.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.model_predict-Tuple{AbstractLinearModelResult, Matrix}","page":"API Reference","title":"Jessamine.model_predict","text":"model_predict(lmr::AbstractLinearModelResult, X::Matrix)\n\nReturn X * coefficients(lmr) + intercept(lmr)\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.model_predict-Tuple{GenomeSpec, Agent{<:Number, <:AbstractGenome, <:AbstractVector, <:AbstractModelResult}, Any}","page":"API Reference","title":"Jessamine.model_predict","text":"model_predict(g_spec::GenomeSpec, agent::Agent, xs; kw_args...)\n\nRun agent.genome on inputs xs and agent.parameter, and form the linear combination of the genome's outputs using the coefficients agent.extra.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.model_symbolic_output-Tuple{Any, Any}","page":"API Reference","title":"Jessamine.model_symbolic_output","text":"model_symbolic_output(g_spec, agent; kw_args...)\n\nBuild a Symbolics form for the output of the final time step of running agent's genome Then use the agent's parameter vector, and feed the symbolic output of the genome as input to model result in agent.extra to make a prediction in symbolic form.\n\nThe kw_args are eventually splatted into model_predict.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.model_symbolic_output-Tuple{GenomeSpec, AbstractGenome, AbstractVector, AbstractModelResult}","page":"API Reference","title":"Jessamine.model_symbolic_output","text":"model_symbolic_output(g_spec, genome, parameter, model_result; kw_args...)\n\nBuild a Symbolics form for the output of the final time step of running genome Then use the parameter vector, the symbolic output of the genome, and feed the symbolic output of the genome as input to model result in agent.extra to make a prediction in symbolic form.\n\nThe kw_args are eventually splatted into model_predict.\n\nReturn a named tuple with lots of useful fields.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.model_symbolic_output-Tuple{NamedTuple}","page":"API Reference","title":"Jessamine.model_symbolic_output","text":"model_symbolic_output(fit_result::NamedTuple; kw_args...)\n\nCall model_symbolic_output with the g_spec::GenomeSpec and best_agent::Agent found in fit_result. Return the result.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.model_sympy_output-Tuple{Any, Any}","page":"API Reference","title":"Jessamine.model_sympy_output","text":"model_sympy_output(g_spec, agent; kw_args...)\n\nBuild a SymPy form for the output of the final time step of running agent's genome Then use the agent's parameter vector, and feed the symbolic output of the genome as input to model result in agent.extra to make a prediction in symbolic form.\n\nKey-word arguments:\n\nassumptions = Dict(:extended_real => true):\n\nAssumptions to use when creating symbols.\n\nOther kw_args are eventually splatted into model_predict.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.model_sympy_output-Tuple{GenomeSpec, AbstractGenome, AbstractVector, AbstractModelResult}","page":"API Reference","title":"Jessamine.model_sympy_output","text":"model_sympy_output(g_spec, genome, parameter, model_result; kw_args...)\n\nBuild a SymPy form for the output of the final time step of running genome Then use the parameter vector, the symbolic output of the genome, and feed the symbolic output of the genome as input to model result in agent.extra to make a prediction in symbolic form.\n\nKey-word arguments:\n\nassumptions = Dict(:extended_real => true):\n\nAssumptions to use when creating symbols.\n\nOther kw_args are eventually splatted into model_predict.\n\nReturn a named tuple with lots of useful fields.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.mutate","page":"API Reference","title":"Jessamine.mutate","text":"mutate(rng::AbstractRNG, m_dist::MutationDist, x)\n\nRandomly change x, using random numbers supplied by RNG and the distributions specified by m_dist.\n\n\n\n\n\n","category":"function"},{"location":"reference/api/#Jessamine.new_genome-Tuple{Random.AbstractRNG, SelectionDist, MutationDist, Population}","page":"API Reference","title":"Jessamine.new_genome","text":"new_genome(rng::AbstractRNG, s_dist::SelectionDist, m_dist::MutationDist, pop::Population)\n\nProduce a new offspring genome. Pick two parents using tournament selection. Recombine their genomes, and apply mutations.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.next_generation-Tuple{Random.AbstractRNG, GenomeSpec, SelectionDist, MutationDist, Population, Any}","page":"API Reference","title":"Jessamine.next_generation","text":"next_generation(\n    rng::AbstractRNG,\n    g_spec::GenomeSpec,\n    s_dist::SelectionDist,\n    m_dist::MutationDist,\n    pop::Population,\n    grow_and_rate;\n    sense = MinSense)\n\nProduce the next generation of a population by selection and mutation.  The offsprings' genomes are produced by new_genome, and fed to grow_and_rate(rng, g_spec, genome), which should \"grow\" each organism and return the rating and extra data as an Agent, which is inserted into the population.  The sense parameter specifies whether selection should minimize or maximize the rating.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.num_instructions-Tuple{Genome}","page":"API Reference","title":"Jessamine.num_instructions","text":"num_instructions(genome::Genome)\n\nReturn the total number of instructions in all blocks in the genome.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.num_operands-Tuple{AbstractVector}","page":"API Reference","title":"Jessamine.num_operands","text":"num_operands(xs::AbstractVector)\n\nReturn the total number of instruction operands in xs.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.num_operands-Tuple{Genome}","page":"API Reference","title":"Jessamine.num_operands","text":"num_operands(genome::Genome)\n\nReturn the total number of operands in all instructions in the genome.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.num_operands-Tuple{Instruction}","page":"API Reference","title":"Jessamine.num_operands","text":"num_operands(instruction::Instruction)\n\nReturn the total number of operands in the instruction.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.op_eval","page":"API Reference","title":"Jessamine.op_eval","text":"op_eval(op::AbstractGeneOp, operands::AbstractVector)\n\nReturn op applied to a list of operands.\n\n\n\n\n\n","category":"function"},{"location":"reference/api/#Jessamine.op_eval_add_into!-Tuple{AbstractVector, AbstractGeneOp, AbstractVector}","page":"API Reference","title":"Jessamine.op_eval_add_into!","text":"op_eval_add_into!(dest, op, operands)\n\nApply the operator op to operands and add the result elementwise in place to dest.  The default is to do dest .= dest .+ op_eval(op, operands).\n\nThe idea here is that in-place accumulation instead of the default implementation can sometimes avoid the allocation of a temporary array to hold the result of op_eval(...).\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.parameter_solver_optimization_function-Tuple{AbstractSolverSpec, Any}","page":"API Reference","title":"Jessamine.parameter_solver_optimization_function","text":"parameter_solver_optimization_function(sol_spec, f)\n\nBuild an OptimizationFunction around the function f(genome_parameter_vector, _). The resulting objective function object will eventually be passed to solve(). The default implementation returns OptimizationFunction(f).\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.parameter_solver_optimization_problem-Tuple{AbstractSolverSpec, GenomeSpec, Any, Any}","page":"API Reference","title":"Jessamine.parameter_solver_optimization_problem","text":"parameter_solver_optimization_problem(sol_spec, g_spec, optim_fn, context)\n\nBuild an OptimizationProblem around the OptimizationFunction. The resulting problem object will be passed to solve. The default implementation returns OptimizationProblem(optim_fn, zeros(...), context, sense=MinSense).\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.parameter_solver_solve-Tuple{AbstractSolverSpec, Any}","page":"API Reference","title":"Jessamine.parameter_solver_solve","text":"parameter_solver_solve(sol_spec, optim_prob)\n\nSolve the problem generated by parameter_solver_optimization_problem(). The default implementation returns solve(optim_prob, NelderMeade())\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.pick_parent-Tuple{Random.AbstractRNG, SelectionDist, Population}","page":"API Reference","title":"Jessamine.pick_parent","text":"pick_parent(rng::AbstractRNG, s_dist::SelectionDist, pop::Population)::Agent\n\nUsing tournament selection, pick a parent from the population.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.prediction_performance-Tuple{AbstractMachineSpec, AbstractVector, AbstractVector}","page":"API Reference","title":"Jessamine.prediction_performance","text":"prediction_performance(mn_spec, y_hat, y_ref)\n\nReturn a performance measure of the prediction y_hat compared to reference y_ref.\n\nThe default implementation applies the callable mn_spec.performance(y_hat, y_ref)\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.random_genome-Tuple{Random.AbstractRNG, GenomeSpec, MutationDist, Distributions.Distribution}","page":"API Reference","title":"Jessamine.random_genome","text":"random_genome(rng::AbstractRNG, g_spec::GenomeSpec, m_spec::MutationSpec, arity_dist::Distribution)\n\nProduce a random genome.  It will contain an instruction block for each mutable slot in the state array as specified by g_spec. Each block will contain one random instruction.  The operator is picked uniformly at random from m_spec.op_inventory.  The number of operands is drawn from arity_dist.  The operands are drawn uniformly from the set of possible indices.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.random_initial_population-Tuple{Random.AbstractRNG, EvolutionSpec, Distributions.Distribution}","page":"API Reference","title":"Jessamine.random_initial_population","text":"random_initial_population(\n    rng::AbstractRNG,\n    e_spec::EvolutionSpec,\n    arity_dist::Distribution;\n    sense = MinSense)\n\nInitialize a random initial population from an e_spec.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.random_initial_population-Tuple{Random.AbstractRNG, GenomeSpec, MutationDist, Distributions.Distribution, SelectionSpec, Any}","page":"API Reference","title":"Jessamine.random_initial_population","text":"random_initial_population(\n    rng::AbstractRNG,\n    g_spec::GenomeSpec,\n    m_dist::MutationDist,\n    arity_dist::Distribution,\n    s_spec::SelectionSpec,\n    grow_and_rate;\n    sense = MinSense)::Population\n\nMake a random initial population.  The number of genomes is specified by adding the number of new genomes per generation and the number of keepers given in s_spec.  Other parameters are passed to random_genome to produce random genomes. The grow_and_rate function is called with rng, g_spec, genome.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.replace_near_integer-Tuple{SymbolicUtils.BasicSymbolic}","page":"API Reference","title":"Jessamine.replace_near_integer","text":"replace_near_integer(expr; tolerance=1.0e-10)\n\nRound literal floating-point numbers within a symbolic expression. Specifically, if a number a differs round(a) by less thant tolerance, it gets replaced by round(a).\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.run_genome-Tuple{GenomeSpec, AbstractGenome, AbstractArray, Any}","page":"API Reference","title":"Jessamine.run_genome","text":"run_genome(g_spec::GenomeSpec, genome::Genome, parameter::AbstractVector, input)::Vector\n\nBuild a work space vector using zeros for each output and scratch slot, followed by the parameter vector, then the input vector.  Evaluate the instructions in genome, repeating the evaluation g_spec.num_time_steps.  Return an array that contains, for each time step, the elements 1 through g_spec.output_size of the work space vector.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.run_genome_symbolic-Tuple{GenomeSpec, AbstractGenome}","page":"API Reference","title":"Jessamine.run_genome_symbolic","text":"run_genome_symbolic(g_spec, genome; paramter_sym=:p, input_sym=:x)\n\nBuild a symbolic form for the output of the final time step of running genome.  The parameter vector and input vector are Symbolics objects of the form p[j] and x[j].  The variable names can be specified with the keyword arguments.\n\nReturn a named tuple (p, x, z) where p and x, are vectors of Symbolics objects used to represent genome parameters and inputs; and w is a vector of genome outputs in symbolic form.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.short_show-Tuple{Any}","page":"API Reference","title":"Jessamine.short_show","text":"short_show([io::IO], x)\n\nPrint a short version of x to io, using stdout by default.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.show_symbolic-Tuple{GenomeSpec, Genome}","page":"API Reference","title":"Jessamine.show_symbolic","text":"show_symbolic(g_spec, genome; output_sym, scratch_sym, parameter_sym, input_sym)\n\nReturn a matrix with three columns. The first is row indices 1, 2, etc. The second is the inital worspace vector in symbolic form. The third is the future workspace state in symbolic form.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.show_symbolic-Tuple{NamedTuple}","page":"API Reference","title":"Jessamine.show_symbolic","text":"show_symbolic(fit_result::NamedTuple; kw_args...)\n\nCall show_symbolic with the g_spec::GenomeSpec and best_agent.genome::AbstractGenome found in fit_result. Return the result.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.splat_or_default-Tuple{Any, Any, Any}","page":"API Reference","title":"Jessamine.splat_or_default","text":"splat_or_default(op, def, operands)\n\nReturn the result of applying op to the operands, with op([]) = def and op([x1...]) = op(x1...). The Base.splat function doesn't have a way to deal with the first case.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.vns_keep_going-Tuple{AbstractPopulationCondition}","page":"API Reference","title":"Jessamine.vns_keep_going","text":"vns_keep_going(::AbstractPopulationCondition)\n\nReturn whether VNS should continue to another epoch. Utility function.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#Jessamine.workspace_size-Tuple{GenomeSpec}","page":"API Reference","title":"Jessamine.workspace_size","text":"workspace_size(g_spec::GenomeSpec)\n\nReturn the number of elements in the workspace vector specified by g_spec.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#StatisticalTraits.supports_weights-Tuple{Type{JessamineDeterministic}}","page":"API Reference","title":"StatisticalTraits.supports_weights","text":"supports_weights(::Type{JessamineModel})\n\nReturn true.  A Jessamine model supports weight vectors for training data if the inner model does.  Since MLJ determines this at the type level, it has to return true here to allow a weight vector to be passed in.  The fit implementation throws an exception if a non-nothing weight vector is passed in and the inner model does not allow a weight vector.\n\n\n\n\n\n","category":"method"},{"location":"reference/api/#StatisticalTraits.supports_weights-Tuple{Type{JessamineProbabilistic}}","page":"API Reference","title":"StatisticalTraits.supports_weights","text":"supports_weights(::Type{JessamineModel})\n\nReturn true.  A Jessamine model supports weight vectors for training data if the inner model does.  Since MLJ determines this at the type level, it has to return true here to allow a weight vector to be passed in.  The fit implementation throws an exception if a non-nothing weight vector is passed in and the inner model does not allow a weight vector.\n\n\n\n\n\n","category":"method"},{"location":"#Jessamine.jl","page":"Jessamine.jl","title":"Jessamine.jl","text":"","category":"section"},{"location":"","page":"Jessamine.jl","title":"Jessamine.jl","text":"Documentation for Jessamine.jl","category":"page"},{"location":"","page":"Jessamine.jl","title":"Jessamine.jl","text":"","category":"page"},{"location":"#Introduction","page":"Jessamine.jl","title":"Introduction","text":"","category":"section"},{"location":"","page":"Jessamine.jl","title":"Jessamine.jl","text":"Note Jessamine is experimental and under development. Expect the source code to change in unpredictable ways.","category":"page"},{"location":"#Index","page":"Jessamine.jl","title":"Index","text":"","category":"section"},{"location":"","page":"Jessamine.jl","title":"Jessamine.jl","text":"","category":"page"}]
}
