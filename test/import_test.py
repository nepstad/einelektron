import sys
import unittest
sys.path.append("..")


class TestImportAllEinPartikkelModules(unittest.TestCase):
	"""
	Test that all EinPartikkel submodules, classes and functions may be imported
	"""
	
	def setUp(self):
		einpartikkel = __import__("einpartikkel")
		self.globals = {"einpartikkel": einpartikkel}
	
	def IsModule(self, obj):
		if hasattr(obj, "__all__"):
			return True
		return False
	
	def GetModuleMembers(self, module):
		assert self.IsModule(module)
		return module.__all__
	
	def ImportMembers(self, module, curIndent=""):
		print "%sNow importing members of %s" % (curIndent, module.__name__)
		for m in self.GetModuleMembers(module):
			
			def importModule(parent, module):
				print "%s  importing %s from %s" % (curIndent, m, parent.__name__)
				importedParentModule = __import__(parent.__name__, globals=self.globals, fromlist=[module])
				return importedParentModule.__getattribute__(module)
			
			self.assert_(importModule(module, m))
			importedModule = importModule(module, m)
			if self.IsModule(importedModule):
				self.globals[m] = importedModule
				self.ImportMembers(importedModule, curIndent+"  ")
			
			
	def test_import_members(self):
		self.ImportMembers(__import__("einpartikkel"))


if __name__ == "__main__":
	unittest.main()
