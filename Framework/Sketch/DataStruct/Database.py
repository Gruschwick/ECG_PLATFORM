

class Database(object):
	def __init__(self, path, loader):
		self.path = path	
		self.loader = loader
		self.Procedures = []

	def addProcedure(self,procedure):
		self.Procedures.append(procedure)

	def load(self,path):
		return self.loader.DataFrom(path)

	def DataFields(self):
		return self.loader.DataFields()
