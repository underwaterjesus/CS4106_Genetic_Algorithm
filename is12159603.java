import java.util.*;
import java.io.*;
import java.lang.Math;
import javax.swing.JFrame;
import javax.swing.JPanel;
import java.awt.Graphics;
import java.awt.GridLayout;
import javax.swing.BoxLayout;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;

public class is12159603
{
	private static final String INPUT = "input.txt";
	private static final int MAX = 100;
	public static int c = 0;
	private static File file;
	private static Scanner scanner;
	
	private static JFrame frame;
    private static JPanel panel;
    private static JTextField jPopulation;
    private static JTextField jGenerations;
	private static JTextField jCrossover;
    private static JTextField jMutation;
	
	private static String sPopulation;
	private static String sGenerations;
	private static String sCrossover;
	private static String sMutation;
	
	private static int N;
	private static int P;
	private static int G;
	private static int Cr;
	private static int Mu;
	
	private static int[][] adjacencyMatrix;
	private static int[][] currentPopulation;
	private static int[][] nextPopulation;
	
	private static GraphVisualization graph;
	
	public static void main(String args[])
	{
		/************************************
		 ** Create input panes, verify input
		 ** and parse input. Give appropriate
		 ** error messages.
		 ************************************/
		try
		{
			panel = new JPanel();
			panel.setLayout(new GridLayout(0, 4, 2, 2));
			
			jPopulation = new JTextField();
			jGenerations = new JTextField();
			jCrossover = new JTextField();
			jMutation = new JTextField();
			
			panel.add(new JLabel("Population(P):"));
			panel.add(jPopulation);
			
			panel.add(new JLabel("Generations(G)"));
			panel.add(jGenerations);
			
			panel.add(new JLabel("Crossover Rate(Cr)"));
			panel.add(jCrossover);
			
			panel.add(new JLabel("Mutation Rate(Mu)"));
			panel.add(jMutation);
			
			int option = JOptionPane.showConfirmDialog(frame, panel, "Genetic details:",
														JOptionPane.YES_NO_OPTION, JOptionPane.INFORMATION_MESSAGE);
														
			if (option == JOptionPane.YES_OPTION)
			{
				sPopulation = jPopulation.getText();
				sGenerations = jGenerations.getText();
				sCrossover =  jCrossover.getText();
				sMutation = jMutation.getText();
				
				if(verifyInt(sPopulation))
				{
					P = Integer.parseInt(sPopulation.trim());
				}
				else
				{
					JOptionPane.showMessageDialog(null, "Error parsing population size", "Error", 1);
					System.exit(0);
				}
				
				if(verifyInt(sGenerations))
				{
					G = Integer.parseInt(sGenerations.trim());
				}
				else
				{
					JOptionPane.showMessageDialog(null, "Error parsing number of generations", "Error", 1);
					System.exit(0);
				}
				
				if(verifyInt(sCrossover))
				{
					Cr = Integer.parseInt(sCrossover.trim());
				}
				else
				{
					JOptionPane.showMessageDialog(null, "Error parsing crossover rate", "Error", 1);
					System.exit(0);
				}
				
				if(verifyInt(sMutation))
				{
					Mu = Integer.parseInt(sMutation.trim());
				}
				else
				{
					JOptionPane.showMessageDialog(null, "Error parsing mutation rate", "Error", 1);
					System.exit(0);
				}
				
				if( (P < 1) || (G < 1) || Cr < 0  || Cr > 100 || Mu < 0 || Mu > 100 || (Cr + Mu) > 100 )
				{
					JOptionPane.showMessageDialog(null, "Forbidden values entered. Please remember:\n" +
															" -P and G must be positive\n" +
															" -Cr and Mu must be between (0 - 100)\n" +
															" -The sum of Cr and Mu must be between (0 - 100)", "Error", 1);
					
					System.exit(0);
				}
			}
			else
			{
				
				System.exit(0);
			}
		}
		catch(Exception e)
		{
			System.out.println(e);
			System.exit(0);
		}
		
		/************************************
		 ** Read input file.
		 ************************************/
		try
		{
			file = new File(INPUT);
			scanner = new Scanner(file);
		}
		catch(Exception e)
		{
			System.out.println(e);
			System.exit(0);
		}
		
		N = 0;
		ArrayList<int[]> list = new ArrayList<int[]>();
		
		/************************************
		 ** Create arraylist of Arrays. Each
		 ** array represents a line of the
		 ** file. This is to avoid reading
		 ** the file twice. 
		 ************************************/
		//try
		//{
			while(scanner.hasNext())
			{
				String[] line = (scanner.nextLine()).split(" ");
				int[] arr = { (Integer.parseInt(line[0])), (Integer.parseInt(line[1])) };
				
				list.add(arr);
				
				if(arr[1] > N)
					N = arr[1];
			}
			
			N++;
			adjacencyMatrix = new int[N][N];
			
			/************************************
			** Populate adjacency matrix from
			** arraylist of arrays.
			************************************/
			for(int[] i : list)
			{
				adjacencyMatrix[ i[0] ][ i[1] ] = 1;
				adjacencyMatrix[ i[1] ][ i[0] ] = 1;
			}
			
			/************************************
			** Print adjacency matrix.
			** Asked for in specification.
			************************************/
			System.out.println(" Adjacency Matrix:");
			for(int i = 0; i < N; i++)
			{
				System.out.print(" |");
				for(int j = 0; j < N; j++)
					System.out.print(" " + ( adjacencyMatrix[i][j]) + " |");
				System.out.println();
			}
			
			/************************************
			** Declare 2-d arrays.
			** Asked for in specification.
			************************************/
			currentPopulation = new int[P][N];
			nextPopulation = new int[P][N];
			
			/************************************
			** Populate G(0) with P unique 
			** random orderings
			************************************/
			for(int i =0; i < P; i++)
			{
				for(int j = 0; j < N; j++)
				{
					currentPopulation[i][j] = j;
				}
				for(int k = 1; k < N; k++)
				{
					int r = (int)(Math.random() * (k + 1));
					int t = currentPopulation[i][k];
					
					currentPopulation[i][k] = currentPopulation[i][r];
					currentPopulation[i][r] = t;
				}
				for(int m = 0; m < i; m++)
				{
					if(Arrays.equals( currentPopulation[i], currentPopulation[m] ))
					{
						i--; break;
					}
				}
			}
			
			/************************************
			** Print orderings in G(0) 
			************************************/
			System.out.println();
			System.out.println(" Initial population of orderings G(0):");
			
			for(int i = 0; i < P; i++)
			{
				System.out.print(" #" + i + "-");
				System.out.println(" " + Arrays.toString(currentPopulation[i]));
			}
			System.out.println();
			
			//graph = new GraphVisualization(adjacencyMatrix, currentPopulation[0], N);
			//graph = new GraphVisualization(adjacencyMatrix, currentPopulation[P - 1], N);
			System.out.println();
			
			/**************************************
			** s is used later to dived population
			** into 3 parts. In specification.
			***************************************/
			int s = (P / 3);
			
			/************************************
			** Loop through generations
			************************************/
			for(int g = 0; g < G; g++)
			{
				/************************************
				** Sort population. Best at index 0.
				************************************/
				sortOrderings( currentPopulation, 0, (P - 1) );
				
				System.out.println("Best performing ordering from G(" + g + "):");
				System.out.println("  " + Arrays.toString(currentPopulation[0]));
				
				/************************************
				** Display best performer of
				** generation
				************************************/
				graph = new GraphVisualization(adjacencyMatrix, currentPopulation[0], N, g);
				
				/************************************
				** Copy best performing third into
				** bottom performing third of
				** population. In specification.
				************************************/
				for(int i = 0; i < s; i++)
					currentPopulation[ (P - s) + i ] = currentPopulation[i].clone();
				
				int p = P;
				
				for(int i = 0; i < P; i++)
				{
					int Pr = (int)(Math.random() * (MAX + 1));
					
					/************************************
					** Perform crossover, else perform
					** peform mutation, else straight
					** copy to next generation.
					************************************/
					if( (Cr >= Pr) && (p > 1) )
					{
						int r = (int)(Math.random() * (p));
						
						while(currentPopulation[r] == null)
						{
							r++;
						}
						
						int[] tempArr = currentPopulation[r];
						currentPopulation[r] = null;
						
						r = (int)(Math.random() * (--p));
						
						while(currentPopulation[r] == null)
						{
							r++;
						}
						
						crossover(tempArr, currentPopulation[r]);
						
						nextPopulation[i] = tempArr;
						nextPopulation[++i] = currentPopulation[r];
						
						currentPopulation[r] = null;
						p--;
					}
					else if( (Cr + Mu) >= Pr )
					{
						int r = (int)(Math.random() * (p));
						
						while(currentPopulation[r] == null)
						{
							r++;
						}
						
						mutate(currentPopulation[r]);
						nextPopulation[i] = currentPopulation[r];
						currentPopulation[r] = null;
						
						p--;
					}
					else
					{
						int r = (int)(Math.random() * (p));
						
						while(currentPopulation[r] == null)
						{
							r++;
						}
						
						nextPopulation[i] = currentPopulation[r];
						currentPopulation[r] = null;
						
						p--;
					}
				}
				
				/************************************
				** Copy new population back into
				** current population.
				************************************/
				for(int i = 0; i < P; i++)
					currentPopulation[i] = nextPopulation[i].clone();
				
				/************************************
				** Test: printing out population
				************************************/
				for(int i = 0; i < P; i++)
				{
				System.out.print(" #" + i + "-");
				System.out.println(" " + Arrays.toString(nextPopulation[i]));
				}
			}
		//}
		//catch(Exception e)
		//{
			//System.out.println(e);
			//System.exit(0);
		//}
	}
	
	public static boolean verifyInt(String test)
	{
		String pattern = "(((\\s))*((-)?)([0-9])+((\\s))*){1}";
		
		return test.matches(pattern);
	}
	
	/************************************
	** Implementation of given fitness
	** cost(total length of lines).
	************************************/
	public static int fitnessCost(int[] ordering)
	{
		double ret = 0.0;
		double radius = 100.0;
		double chunk = ( (Math.PI * 2.0) / ( (double)N ) );
		
		/************************************
		** Loop through upper half of
		** adjacencyMatrix. If edge found
		** work out distance using line
		** functions
		** https://orion.math.iastate.edu/dept/links/formulas/form2.pdf
		************************************/
		for(int i = 0; i < N; i++)
				for(int j = i + 1; j < N; j++)
					if( adjacencyMatrix[ ordering[i] ][ ordering[j] ] == 1 )
					{
						double x1 = ( Math.cos(i * chunk) * radius );
						double y1 = ( Math.sin(i * chunk) * radius );
						double x2 = ( Math.cos(j * chunk) * radius );
						double y2 = ( Math.sin(j * chunk) * radius );
						
						ret += ( Math.sqrt( ( Math.pow( (x2 - x1), 2 ) ) + ( Math.pow( (y2 - y1), 2 ) ) ) );
					}
		
		return (int)ret;
	}
	
	/************************************
	** Custom mergesort implementation
	************************************/
	public static void sortOrderings(int[][] population, int low, int hi)
	{
		if (low >= hi)
			return;
		
		int mid = (low + hi) / 2;
		
		sortOrderings(population, low, mid);
		sortOrderings(population, (mid + 1), hi);
		
		mergeParts(population, low, mid, hi);
	}
	
	public static void mergeParts(int[][]part, int low, int mid, int hi)
	{
		int left = low; int right = mid + 1; int temp = 0;
		int numElements = (hi - low) + 1;
		int[][] tempArr = new int[part.length][N];
		
		while(left <= mid && right <= hi)
			if( fitnessCost(part[left]) <= fitnessCost(part[right]) )
				tempArr[temp++] = part[left++];
			else
				tempArr[temp++] = part[right++];
		
		while(left <= mid)
		tempArr[temp++] = part[left++];
		while(right <= hi)
		tempArr[temp++] = part[right++];
		
		for(int i = 0 ; i < numElements ; i++)
			part[low + i] = tempArr[i];
	
	}
	
	/************************************
	** Mutate population member. Pick 2
	** random indices and swap them.
	************************************/
	public static void mutate(int[] gene)
	{
		int r = (int)(Math.random() * (gene.length));
		int r2 = (int)(Math.random() * ( (gene.length) - 1));
		
		while(r2 == r)
			r2++;
		
		int t = gene[r];
		gene[r] = gene[r2];
		gene[r2] = t;
	}

	/************************************
	** Crossover two population members.
	************************************/
	public static void crossover(int[] gene1, int[] gene2)
	{
		/************************************
		** Pick random crossover point
		************************************/
		int r = ( (int)(Math.random() * ((gene1.length) - 2)) ) + 1;
		
		/************************************
		** Create two tempory "genes"
		************************************/
		int[] t = new int[gene1.length]; int[] t2 = new int[gene1.length];
		
		/************************************
		** Copy section 1 of gene1 into t2,
		** and section 2 of gene2 into t2.
		************************************/
		System.arraycopy(gene1, 0, t2, 0, r);
		System.arraycopy(gene2, r, t2, r, ((gene1.length) - r));
		
		/************************************
		** Copy section 1 of gene2 into t,
		** and section 2 of gene1 into t.
		************************************/
		System.arraycopy(gene2, 0, t, 0, r);
		System.arraycopy(gene1, r, t, r, ((gene1.length) - r));
		
		/************************************
		** Copy new "genes" back into
		** originals
		************************************/
		System.arraycopy(t, 0, gene1, 0, gene1.length);
		System.arraycopy(t2, 0, gene2, 0, gene1.length);
		
		/************************************
		** Rest of method removes duplicates
		** from gene1 and gene2
		************************************/
		boolean[] dupes = new boolean[gene1.length];
		
		for(int i = 0 ; i < r; i++)
			for( int j = 0 ; j < r; j++)
			{
				if(gene1[i] == gene2[j])
				{
					dupes[ gene1[i] ] = true;
					break;
				}
			}
		
		int i = 0; int j = 0;
		while( i < r && j < r)
		{
			if(!( dupes[ gene1[i] ] ))
			{
				while( dupes[ gene2[j] ] )
					j++;
				
				int tmp = gene1[i];
				gene1[i] = gene2[j];
				gene2[j] = tmp;
				j++;
			}
			
			i++;
		}
	}
	
	/************************************
	** Line drawing method given in
	** specification.
	************************************/
	public static class GraphVisualization extends JFrame
	{
		private static final String TITLE = "Graph Visualization";
		private static final int WIDTH = 960;
		private static final int HEIGHT = 960;
		
		private int[][] adjacencyMatrix;
		private int[] ordering;
		
		private int numberOfVertices;
		private double chunk;
		
		public GraphVisualization(int[][] adjacencyMatrix, int[] ordering, int numberOfVertices, int gen)
		{
			this.adjacencyMatrix = adjacencyMatrix;
			this.ordering = ordering;
			this.numberOfVertices = numberOfVertices;
			this.chunk = ( (Math.PI * 2) / ( (double)numberOfVertices ) );
			
			setTitle(TITLE + ": G(" + gen + ")");
			setSize(WIDTH, HEIGHT);
			setVisible(true);
			setDefaultCloseOperation(EXIT_ON_CLOSE);
		}
		
		@Override
		public void paint(Graphics g)
		{

			
			int radius = 100;
			int mov = 200;
			
			for(int i = 0; i < numberOfVertices; i++)
				for(int j = i + 1; j < numberOfVertices; j++)
					if( adjacencyMatrix[ ordering[i] ][ ordering[j] ] == 1 )
					{
						g.drawLine( (int)( Math.cos(i * chunk) * radius ) + mov,
									(int)( Math.sin(i * chunk) * radius ) + mov,
									(int)( Math.cos(j * chunk) * radius ) + mov,
									(int)( Math.sin(j * chunk) * radius ) + mov);
					}
		}
	}
}