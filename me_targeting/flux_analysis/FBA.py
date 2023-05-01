'''
Created on 2014. 6. 4.

@author: user
'''

from gurobipy import *

# import Simulator
from me_targeting.flux_analysis import Simulator


class FBA(Simulator.Simulator):
	'''
	classdocs
	'''

	def __init__(self):
		'''
		Constructor
		'''
		super(FBA, self).__init__()

	def run_FBA(self, new_objective='', flux_constraints={}, inf_flag=False, internal_min_flag=False, mode='max'):
		# print 'Start simple FBA simulation ... '

		model_metabolites = self.model_metabolites
		model_reactions = self.model_reactions

		if new_objective == '':
			objective = self.objective
		else:
			objective = new_objective

		Smatrix = self.Smatrix

		LowerFlux = self.Lower_Boundary_Constraint
		UpperFlux = self.Upper_Boundary_Constraint

		if inf_flag == False:
			for key in LowerFlux.keys():
				if LowerFlux[key] == float("-inf"):
					LowerFlux[key] = -1000.0

			for key in UpperFlux.keys():
				if UpperFlux[key] == float("inf"):
					UpperFlux[key] = 1000.0

		pairs, coffvalue = multidict(Smatrix)
		pairs = tuplelist(pairs)

		m = Model('FBA')
		m.setParam('OutputFlag', 0)
		m.reset()

		# create variables
		v = {}
		fplus = {}
		fminus = {}


		m.update()

		for each_reaction in model_reactions:
			if each_reaction in flux_constraints.keys():
				v[each_reaction] = m.addVar(lb=flux_constraints[each_reaction][0],
				                            ub=flux_constraints[each_reaction][1], name=each_reaction)
			else:
				v[each_reaction] = m.addVar(lb=LowerFlux[each_reaction], ub=UpperFlux[each_reaction],
				                            name=each_reaction)
			fplus[each_reaction] = m.addVar(lb=0.0, ub=1000.0, name=each_reaction)
			fminus[each_reaction] = m.addVar(lb=0.0, ub=1000.0, name=each_reaction)
		m.update()

		# Add constraints
		for each_metabolite in model_metabolites:
			if len(pairs.select(each_metabolite, '*')) == 0:
				continue
			m.addConstr(quicksum(v[reaction] * coffvalue[metabolite, reaction] for metabolite, reaction in
			                     pairs.select(each_metabolite, '*')) == 0)

		m.update()
		if mode == 'max':
			m.setObjective(v[objective], GRB.MAXIMIZE)
		elif mode == 'min':
			m.setObjective(v[objective], GRB.MINIMIZE)

		m.optimize()
		print(m.status)
		if m.status == 2:
			obj_value = m.ObjVal
			print('Objective value : %s'%m.ObjVal)
			if internal_min_flag == True:
				m.addConstr(fplus[objective]-fminus[objective] == obj_value)

				m.addConstr(quicksum((fplus[reaction]-fminus[reaction]) * coffvalue[metabolite, reaction] for metabolite, reaction in
			                     pairs.select(each_metabolite, '*')) == 0)

				for each_reaction in model_reactions:
					m.addConstr(fplus[each_reaction]-fminus[each_reaction] == v[each_reaction])

				m.update()
				m.setObjective(quicksum((fplus[each_reaction]+fminus[each_reaction]) for each_reaction in model_reactions),
				               GRB.MINIMIZE)
				m.optimize()
				if m.status == 2:
					obj_value = m.ObjVal
					print('Flux minimization objective value : %s'%m.ObjVal)
					ReactionFlux = {}
					for reaction in model_reactions:
						ReactionFlux[reaction] = float(v[reaction].x)
						if abs(float(v[reaction].x)) <= 1e-6:
							ReactionFlux[reaction] = 0.0
					return m.status, obj_value, ReactionFlux
			else:
				ReactionFlux = {}
				for reaction in model_reactions:
					ReactionFlux[reaction] = float(v[reaction].x)
					if abs(float(v[reaction].x)) <= 1e-6:
						ReactionFlux[reaction] = 0.0

				return m.status, obj_value, ReactionFlux
		return m.status, False, False

	def run_pFBA(self, new_objective='', flux_constraints={}, inf_flag=False, mode='max'):
		internal_min_flag=True
		# print 'Start simple FBA simulation ... '

		model_metabolites = self.model_metabolites
		model_reactions = self.model_reactions

		if new_objective == '':
			objective = self.objective
		else:
			objective = new_objective

		Smatrix = self.Smatrix

		LowerFlux = self.Lower_Boundary_Constraint
		UpperFlux = self.Upper_Boundary_Constraint

		if inf_flag == False:
			for key in LowerFlux.keys():
				if LowerFlux[key] == float("-inf"):
					LowerFlux[key] = -1000.0

			for key in UpperFlux.keys():
				if UpperFlux[key] == float("inf"):
					UpperFlux[key] = 1000.0

		pairs, coffvalue = multidict(Smatrix)
		pairs = tuplelist(pairs)

		m = Model('FBA')
		m.setParam('OutputFlag', 0)
		m.reset()

		# create variables
		v = {}
		fplus = {}
		fminus = {}


		m.update()

		for each_reaction in model_reactions:
			if each_reaction in flux_constraints.keys():
				v[each_reaction] = m.addVar(lb=flux_constraints[each_reaction][0],
				                            ub=flux_constraints[each_reaction][1], name=each_reaction)
			else:
				v[each_reaction] = m.addVar(lb=LowerFlux[each_reaction], ub=UpperFlux[each_reaction],
				                            name=each_reaction)
			fplus[each_reaction] = m.addVar(lb=0.0, ub=1000.0, name=each_reaction)
			fminus[each_reaction] = m.addVar(lb=0.0, ub=1000.0, name=each_reaction)
		m.update()

		# Add constraints
		for each_metabolite in model_metabolites:
			if len(pairs.select(each_metabolite, '*')) == 0:
				continue
			m.addConstr(quicksum(v[reaction] * coffvalue[metabolite, reaction] for metabolite, reaction in
			                     pairs.select(each_metabolite, '*')) == 0)

		m.update()
		if mode == 'max':
			m.setObjective(v[objective], GRB.MAXIMIZE)
		elif mode == 'min':
			m.setObjective(v[objective], GRB.MINIMIZE)

		m.optimize()
		print(m.status)
		if m.status == 2:
			obj_value = m.ObjVal
			print('Objective value : %s'%m.ObjVal)
			if internal_min_flag == True:
				m.addConstr(fplus[objective]-fminus[objective] == obj_value)

				m.addConstr(quicksum((fplus[reaction]-fminus[reaction]) * coffvalue[metabolite, reaction] for metabolite, reaction in
			                     pairs.select(each_metabolite, '*')) == 0)

				for each_reaction in model_reactions:
					m.addConstr(fplus[each_reaction]-fminus[each_reaction] == v[each_reaction])

				m.update()
				m.setObjective(quicksum((fplus[each_reaction]+fminus[each_reaction]) for each_reaction in model_reactions),
				               GRB.MINIMIZE)

				m.optimize()
				if m.status == 2:
					obj_value = m.ObjVal
					print('Flux minimization objective value : %s'%m.ObjVal)
					ReactionFlux = {}
					for reaction in model_reactions:
						ReactionFlux[reaction] = float(v[reaction].x)
						if abs(float(v[reaction].x)) <= 1e-6:
							ReactionFlux[reaction] = 0.0
					return m.status, obj_value, ReactionFlux
			else:
				ReactionFlux = {}
				for reaction in model_reactions:
					ReactionFlux[reaction] = float(v[reaction].x)
					if abs(float(v[reaction].x)) <= 1e-6:
						ReactionFlux[reaction] = 0.0

				return m.status, obj_value, ReactionFlux
		return m.status, False, False

if __name__ == '__main__':
	obj = FBA()
	obj.readModel('./iJO1366.xml', False)
	status, ObjVal, ReactionFlux = obj.run_pFBA()
	print(ObjVal)
	for each_reaction in ReactionFlux:
		print(each_reaction,'\t', ReactionFlux[each_reaction])

	# print 'Test'
