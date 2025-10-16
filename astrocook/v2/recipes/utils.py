from typing import Dict, Any, Callable
import importlib

def get_recipe_schema(category: str, recipe_name: str) -> Dict[str, Any]:
    """
    Dynamically loads the schemas from the corresponding recipe module.
    :param category: The recipe category (e.g., 'edit', 'fit').
    :param recipe_name: The name of the recipe method (e.g., 'x_convert').
    """
    try:
        # Construct the full path to the module (e.g., 'astrocook.v2.recipes.edit')
        module_path = f"astrocook.v2.recipes.{category}" 
        recipe_module = importlib.import_module(module_path)
        
        # Get the global schemas dictionary from the module (e.g., EDIT_RECIPES_SCHEMAS)
        schemas_dict = getattr(recipe_module, f"{category.upper()}_RECIPES_SCHEMAS")
        
        # Retrieve the specific recipe schema
        schema = schemas_dict[recipe_name]
        
    except (ImportError, AttributeError, KeyError) as e:
        raise NotImplementedError(f"Schema not found for recipe {category}.{recipe_name}: {e}")
        
    return schema