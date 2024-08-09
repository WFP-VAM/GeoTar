def replace_local_to_s3(data):
    """
  Replaces local file paths with S3 paths in a dictionary.

  Args:
      data: The dictionary containing data with local file paths.

  Returns:
      A new dictionary with S3 paths replacing the local paths.
  """
    new_data = {}
    for key, value in data.items():
        new_value = {}
        for inner_key, inner_value in value.items():
            if inner_key == "mask_shp":
                # Replace local path with S3 path
                base_path = "C:/"
                if inner_value.startswith(base_path):
                    relative_path = inner_value[len(base_path):].replace("\\", "/")
                    new_value[inner_key] = f"s3://geotar.s3.hq/{relative_path}"
                else:
                    new_value[inner_key] = inner_value
            else:
                new_value[inner_key] = inner_value
        new_data[key] = new_value
    return new_data