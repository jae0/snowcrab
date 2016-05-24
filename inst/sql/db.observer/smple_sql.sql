-- The command below will be the simplest to get the object name and type
-- The list will include tables, views and synonyms
-- The command is not case sensitive 

SELECT * FROM cat;

-- The command below will list the name of tables

SELECT table_name FROM user_tables;


-- The command below will list the name of views

SELECT view_name FROM user_views;


-- The command below will list the synonyms and some other info on synonyms

SELECT * FROM user_synonyms;

-- If you want to see the columns in a table, view or synonym
-- The command below will list the columns of the given table/view/synonym
-- replace the name_of_object with the actual name of the table,view or synonym

desc name_of_object;


-- The commands above will give information on the objects
-- To view the actual data from a table, view or synonym
-- Use the command like the one below

SELECT * FROM scs_main;

SELECT tow_date, trawl_number
FROM scs_main
WHERE sex_cd = 0;

