from importlib import import_module

_EXPORTS = {
	"CurveSimulator": (".curvesim", "CurveSimulator"),
	"CurveSimAnimation": (".cs_animation", "CurveSimAnimation"),
	"CurveSimBodies": (".cs_bodies", "CurveSimBodies"),
	"CurveSimBody": (".cs_body", "CurveSimBody"),
	# "CurveSimFluxData": (".cs_flux_data", "CurveSimFluxData"),
	# "csv2df": (".cs_flux_data", "csv2df"),
	"CurveSimLightcurve": (".cs_lightcurve", "CurveSimLightcurve"),
	"CurveSimMCMC": (".cs_mcmc", "CurveSimMCMC"),
	"CurveSimParameters": (".cs_parameters", "CurveSimParameters"),
	"FittingParameter": (".cs_parameters", "FittingParameter"),
	"CurveSimPhysics": (".cs_physics", "CurveSimPhysics"),
	"CurveSimRebound": (".cs_rebound", "CurveSimRebound"),
	"CurveSimResults": (".cs_results", "CurveSimResults"),
}

__all__ = list(_EXPORTS)


def __getattr__(name):
	try:
		module_name, attribute_name = _EXPORTS[name]
	except KeyError as exc:
		raise AttributeError(f"module {__name__!r} has no attribute {name!r}") from exc
	value = getattr(import_module(module_name, __name__), attribute_name)
	globals()[name] = value
	return value

# Old version:
# # __init__.py
# from .curvesim import CurveSimulator
# from .cs_animation import CurveSimAnimation
# from .cs_bodies import CurveSimBodies
# from .cs_body import CurveSimBody
# from .cs_flux_data import CurveSimFluxData, csv2df
# from .cs_lightcurve import CurveSimLightcurve
# from .cs_mcmc import CurveSimMCMC
# from .cs_parameters import CurveSimParameters, FittingParameter
# from .cs_physics import CurveSimPhysics
# from .cs_rebound import CurveSimRebound
# from .cs_results import CurveSimResults
