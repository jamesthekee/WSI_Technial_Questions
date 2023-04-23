select COUNT(DISTINCT(dept_id))
from Project
INNER JOIN Employee 
on Employee.Project_id == Project.Project_id
where Project_name == "maths" -- Substitute maths for the project_name in question