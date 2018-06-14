/*
 *
 * Simple class to convert a MTLR format file to CSV
 *
 */
import java.io.*;

class ConvertDataFiles
{
	public static void main(String args[])
	{
		if (args.length < 3)
			System.out.println("Usage: java ConvertDataFiles convert2MTLR/convert2CSV <infile> <outfile>");

		else
		{
			String infile = args[1];
			String outfile = args[2];

			if (args[0].equals("convert2MTLR"))
			{
				String[][] data = readCSV(infile);
				if (args.length>3)
				{
					if (args[3].equals("FlipCensoredBit"))
					    data = flipCensoredBit(data);
				}
				writeMTLR(data, outfile);
			}
			else if (args[0].equals("convert2CSV"))
			{
				String[][] data = readMTLR(infile);
				if (args.length>3)
				{
					if (args[3].equals("FlipCensoredBit"))
					    data = flipCensoredBit(data);
				}
				writeCSV(data, outfile);
			}
			else
				System.out.println("Usage: java ConvertDataFiles convert2MTLR/convert2CSV <infile> <outfile>");
		}

	}


    public ConvertDataFiles()
    {}


    public static String[][] flipCensoredBit(String[][] raw_data)
    {
		int numRows = raw_data.length;
		int numCols = raw_data[0].length;
		int censorBitIndex=1;

		String[][] data_copy = new String[numRows][numCols];

		for (int i=0; i<numRows; i++)
		{
			for (int j=0; j<numCols; j++)
			{
				if (j!=censorBitIndex)
				    data_copy[i][j] = raw_data[i][j];
				else
				{
					if (raw_data[i][j].equals("1"))
					    data_copy[i][j] = "0";
					else
					    data_copy[i][j] = "1";
				}
			}
		}

		return data_copy;
	}


	public static String[][] readMTLR(String fName)
	{
		String[][] result = null;
		int i, j;

		try
		{
			FileInputStream fstream = new FileInputStream(fName);
			BufferedReader in = new BufferedReader(new InputStreamReader(fstream));

			// read one line to determine how many features there are
			String line = in.readLine();
			String l = line;
			int numFeats=1;
			while (l.contains(" ")) {
				l = l.substring(0, l.lastIndexOf(" ")).trim();
				numFeats++;
			}

            int numLines=1;
			// read all lines to determine how many
			while (in.ready())
			{
				in.readLine();
				numLines++;
			}
			in.close();

			result = new String[numLines][numFeats];

			// reread file to parse data into result array
			fstream = new FileInputStream(fName);
			in = new BufferedReader(new InputStreamReader(fstream));

			for (i=0; i<numLines; i++)
			{
				line = in.readLine();

				for (j=0; j<numFeats-1; j++)
				{
					l = line.substring(0, line.indexOf(" ")).trim();
					line = line.substring(line.indexOf(" ")+1).trim();

					if (l.contains(":"))
					    l = l.substring(l.indexOf(":")+1);

//                    if (l.length() > 0)
						result[i][j] = l;
//					else
//						result[i][j]=" ";
				}
				if (line.contains(":"))
					line = line.substring(line.indexOf(":")+1);
				result[i][j] = line;
			}
			in.close();
		}
		catch (Exception e)
		{
			System.err.println("File input error");
			e.printStackTrace();
		}
		return result;
	}


	public static String[][] readCSV(String fName)
	{
		String[][] result = null;
		int i, j;

		try
		{
			FileInputStream fstream = new FileInputStream(fName);
			BufferedReader in = new BufferedReader(new InputStreamReader(fstream));

			// read one line to determine how many features there are
			String line = in.readLine();
			String l = line;
			int numFeats=1;
			while (l.contains(",")) {
				l = l.substring(0, l.lastIndexOf(",")).trim();
				numFeats++;
			}

            int numLines=0; // first line is header line, so don't count it
			// read all lines to determine how many
			while (in.ready())
			{
				in.readLine();
				numLines++;
			}
			in.close();

			result = new String[numLines][numFeats];

			// reread file to parse data into result array
			fstream = new FileInputStream(fName);
			in = new BufferedReader(new InputStreamReader(fstream));

			// read header line
			in.readLine();

			for (i=0; i<numLines; i++)
			{
				line = in.readLine();

				for (j=0; j<numFeats-1; j++)
				{
					l = line.substring(0, line.indexOf(",")).trim();
					line = line.substring(line.indexOf(",")+1).trim();

					result[i][j] = l;
				}
				result[i][j] = line;
			}
			in.close();
		}
		catch (Exception e)
		{
			System.err.println("File input error");
			e.printStackTrace();
		}
		return result;
	}


	public static void writeMTLR(String[][] data, String fName)
	{
		if (data==null) {
			System.out.println("Failure reading data from input file, cannot write output.");
			return;
		}
		int i, j;

		try
		{
			FileOutputStream out = new FileOutputStream(fName);
			PrintStream p = new PrintStream(out);
			String s;

			// no header
			for (i=0; i<data.length; i++)
			{
				s = "" + data[i][0] + " " + data[i][1];
				for (j=2; j<data[0].length; j++)
					s = s + " " + (j-1) + ":" + data[i][j];
				p.println(s);
			}
			p.close();
		}
		catch (Exception e)
		{
			System.err.println("File output error");
			e.printStackTrace();
		}
	}


	public static void writeCSV(String[][] data, String fName)
	{
		if (data==null) {
			System.out.println("Failure reading data from input file, cannot write output.");
			return;
		}
		int i, j;

		try
		{
			FileOutputStream out = new FileOutputStream(fName);
			PrintStream p = new PrintStream(out);

			// print header
			String s = "time,delta";
			for (j=2; j<data[0].length; j++)
				s = s + ",v" + (j-1);
			p.println(s);

			for (i=0; i<data.length; i++)
			{
				s = "";
				for (j=0; j<data[0].length-1; j++)
					s = s + data[i][j] + ",";
				s = s + data[i][j];
				p.println(s);
			}
			p.close();
		}
		catch (Exception e)
		{
			System.err.println("File output error");
			e.printStackTrace();
		}
	}

}

