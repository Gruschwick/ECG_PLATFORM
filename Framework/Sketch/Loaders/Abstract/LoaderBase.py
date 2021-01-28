
class Base(object):
	
	#return list of attributes loaded from annotations, called by each loader in Comparator
	def DataFields(self):
		pass

	#reutrn BasicData object containing all data possible
	def DataFrom(self,path):
		pass

	pass
