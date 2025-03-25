function updateJson(field, data, jsonPath)

jsonData = fileread(jsonPath);
parsedData = jsondecode(jsonData);
parsedData.(field) = data;
jsonData = jsonencode(parsedData);
imwrite(jsonData, jsonPath)

end