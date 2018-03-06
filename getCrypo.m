function t = getCrypo(filename, dateformat)
    table = readtable(filename);
    table.Var1 = datetime(table.Var1,'InputFormat',dateformat);
    table.Properties.VariableNames(1) = {'Date'};
    t = table;
end
