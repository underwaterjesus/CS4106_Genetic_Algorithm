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
import javax.swing.ButtonGroup;
import javax.swing.JRadioButton;
import javax.swing.SwingUtilities;

public class Test
{
	private static final String INPUT = "input.txt";
	private static final int MAX = 100;
	private static File file;
	private static Scanner scanner;
	
	private static JFrame frame;
    private static JPanel panel;
    private static JTextField jPopulation;
    private static JTextField jGenerations;
	private static JTextField jCrossover;
    private static JTextField jMutation;
    private static JTextField jIterations;
    private static ButtonGroup buttonGroup;
    private static JRadioButton defaultButton;
    private static JRadioButton angGAButton;
	
	private static String sPopulation;
	private static String sGenerations;
	private static String sCrossover;
	private static String sMutation;
	
	private static int I;
	private static int N;
	private static int P;
	private static int G;
	private static int Cr;
	private static int Mu;
	private static int average;
	private static int totalCrossings = 0;
	private static int lowest = Integer.MAX_VALUE;
	private static int highest = 0;
	
	private static double[] betweenessCentralities;
	
	private static int[][] adjacencyMatrix;
	private static int[][] currentPopulation;
	private static int[][] nextPopulation;
	
	private static GraphVisualization graph;

	private static boolean defaultMethod = true;

	private static ArrayList<Integer> listForMedian = new ArrayList<Integer>();
	
	public static void main(String args[])
	{
		try
		{
			panel = new JPanel();
			panel.setLayout(new GridLayout(0, 2, 2, 2));
			
			jPopulation = new JTextField();
			jGenerations = new JTextField();
			jCrossover = new JTextField();
			jMutation = new JTextField();
			jIterations = new JTextField();
			buttonGroup = new ButtonGroup();
			defaultButton = new JRadioButton("Use default fitness function(total edge lengths)", true);
			angGAButton = new JRadioButton("Use AngGA fitness function", false);

			buttonGroup.add(defaultButton);
			buttonGroup.add(angGAButton);
			
			panel.add(new JLabel("Population(P):"));
			panel.add(jPopulation);
			
			panel.add(new JLabel("Generations(G):"));
			panel.add(jGenerations);
			
			panel.add(new JLabel("Crossover Rate(Cr):"));
			panel.add(jCrossover);
			
			panel.add(new JLabel("Mutation Rate(Mu):"));
			panel.add(jMutation);

			panel.add(new JLabel("Iterations:"));
			panel.add(jIterations);

			panel.add(defaultButton);
			panel.add(angGAButton);
			
			int option = JOptionPane.showConfirmDialog(frame, panel, "Genetic details:",
														JOptionPane.YES_NO_OPTION, JOptionPane.INFORMATION_MESSAGE);
														
			if (option == JOptionPane.YES_OPTION)
			{
				sPopulation = jPopulation.getText();
				sGenerations = jGenerations.getText();
				sCrossover =  jCrossover.getText();
				sMutation = jMutation.getText();
				String sIterations = jIterations.getText(); 

				if(defaultButton.isSelected())
					defaultMethod = true;
				else if(angGAButton.isSelected())
					defaultMethod = false;
				
				if(verifyInt(sIterations))
				{
					I = Integer.parseInt(sIterations.trim());
				}
				else
				{
					JOptionPane.showMessageDialog(null, "Error parsing number of iterations", "Error", 1);
					System.exit(0);
				}

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
				
				if( (I < 1) || (P < 1) || (G < 1) || Cr < 0  || Cr > 100 || Mu < 0 || Mu > 100 || (Cr + Mu) > 100 )
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

		try
		{
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
			betweenessCentralities = new double[N];

			for(int[] i : list)
			{
				adjacencyMatrix[ i[0] ][ i[1] ] = 1;
				adjacencyMatrix[ i[1] ][ i[0] ] = 1;
			}
			
			System.out.println(" Adjacency Matrix:");
			for(int i = 0; i < N; i++)
			{
				System.out.print(" |");
				for(int j = 0; j < N; j++)
					System.out.print(" " + ( adjacencyMatrix[i][j]) + " |");
				System.out.println();
			}
			
			if(!defaultMethod)
				getBetweenessCentralities();

			currentPopulation = new int[P][N];
			nextPopulation = new int[P][N];
			
			for(int x = 0; x < I; x++)
			{
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

				//System.out.println();
				//System.out.println(" Initial population of orderings G(0):");
				
				//for(int i = 0; i < P; i++)
				//{
				//	System.out.print(" #" + i + "-");
				//	System.out.println(" " + Arrays.toString(currentPopulation[i]));
				//}
				//System.out.println();
				
				//System.out.println();

				int s = (P / 3);

			
				for(int g = 0; g < G; g++)
				{
					sortOrderings( currentPopulation, 0, (P - 1) );
					
					//System.out.println("Best performing ordering from G(" + g + "):");
					//System.out.println("  " + Arrays.toString(currentPopulation[0]) + "\n");

					for(int i = 0; i < s; i++)
						currentPopulation[ (P - s) + i ] = currentPopulation[i].clone();
					
					int p = P;
					
					for(int i = 0; i < P; i++)
					{
						int Pr = (int)(Math.random() * (MAX + 1));

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
					
					for(int i = 0; i < P; i++)
						currentPopulation[i] = nextPopulation[i].clone();
				}

				sortOrderings( currentPopulation, 0, P - 1 );
				int add = countCrossings(currentPopulation[0]);
				totalCrossings += add;
				highest = Math.max(highest, add);
				lowest = Math.min(lowest, add);
				listForMedian.add(add);
				//graph = new GraphVisualization(adjacencyMatrix, currentPopulation[0], N, G);
			}
			
			//sortOrderings( currentPopulation, 0, P - 1 );
			//int add = countCrossings(currentPopulation[0]);
			//totalCrossings += add;
			//highest = Math.max(highest, add);
			//lowest = Math.min(lowest, add);
			//System.out.println("Best performing ordering from final generation:");
			//System.out.println("  " + Arrays.toString(currentPopulation[0]));
			graph = new GraphVisualization(adjacencyMatrix, currentPopulation[0], N, G);
			//System.out.println("\nCrossings in final ordering: " + countCrossings(currentPopulation[0]));
			Collections.sort(listForMedian); int halfway = (listForMedian.size() / 2) + 1;
			double median = listForMedian.size() % 2 == 0 ? ((double) listForMedian.get(halfway) + (double) listForMedian.get(halfway - 1)) / 2.0 : (double) listForMedian.get(halfway);
			System.out.printf("\nTests:\t\t\t%d\nMean# Crossings:\t%f\nMedian# Crossings:\t%f\nFewest Crossings:\t%d\nMost Crossings:\t\t%d", I, (double)((double)totalCrossings/(double)I), median, lowest, highest);		
		}
		catch(Exception e)
		{
			System.out.println(e.getMessage());
			System.exit(0);
		}
	}
	
	private static boolean verifyInt(String test)
	{
		String pattern = "(((\\s))*((-)?)([0-9])+((\\s))*){1}";
		
		return test.matches(pattern);
	}
	
	private static double fitnessCost(int[] ordering)
	{
		if(!defaultMethod)
		{
			double ret = 0.0;
			double radius = 100.0;
			double chunk = ( (Math.PI * 2.0) / ( (double)N ) );
			int mov = 200;
			
			for(int i = 0; i < N; i++)
			{
				for(int j = 0; j < N; j++)
				{
					if( adjacencyMatrix[ ordering[i] ][ ordering[j] ] == 1 )
					{
						for(int k = j + 1; k < N; k++)
						{
							if( adjacencyMatrix[ ordering[i] ][ ordering[k] ] == 1 )
							{
								double vi_x = ( Math.cos(i * chunk) * radius ) + mov;
								double vi_y = ( Math.sin(i * chunk) * radius ) + mov;
								double vj_x = ( Math.cos(j * chunk) * radius ) + mov;
								double vj_y = ( Math.sin(j * chunk) * radius ) + mov;
								double vk_x = ( Math.cos(k * chunk) * radius ) + mov;
								double vk_y = ( Math.sin(k * chunk) * radius ) + mov;
								
								double m_ij = (vj_y - vi_y) / (vj_x - vi_x); if(m_ij == Double.POSITIVE_INFINITY || m_ij == Double.NEGATIVE_INFINITY)m_ij = 0.0;
								double m_ik = (vk_y - vi_y) / (vk_x - vi_x); if(m_ik == Double.POSITIVE_INFINITY || m_ik == Double.NEGATIVE_INFINITY)m_ik = 0.0;
								
								double d_ij = ( Math.sqrt( ( Math.pow( (vj_x - vi_x), 2 ) ) + ( Math.pow( (vj_y - vi_y), 2 ) ) ) );
								double d_ik = ( Math.sqrt( ( Math.pow( (vk_x - vi_x), 2 ) ) + ( Math.pow( (vk_y - vi_y), 2 ) ) ) );

								double angle = Math.abs( Math.toDegrees( Math.atan( ( m_ij - m_ik ) / ( 1.0 + ( m_ij * m_ik ) ) ) ) );
								if(angle > 180.0)
									angle = 360.0 - angle;
								
								ret += ( angle * d_ij * d_ik * betweenessCentralities[i] * betweenessCentralities[j] * betweenessCentralities[k] );
							}
						}
					}
				}
			}
			
			return ret;
		}
		else
		{
			double ret = 0.0;
			double radius = 100.0;
			double chunk = ( (Math.PI * 2.0) / ( (double)N ) );

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
			
			return ret;
		}
	}
	
	private static void sortOrderings(int[][] population, int low, int hi)
	{
		if (low >= hi)
			return;
		
		int mid = (low + hi) / 2;
		
		sortOrderings(population, low, mid);
		sortOrderings(population, (mid + 1), hi);
		
		mergeParts(population, low, mid, hi);
	}
	
	private static void mergeParts(int[][]part, int low, int mid, int hi)
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
	
	private static void mutate(int[] gene)
	{
		int r = (int)(Math.random() * (gene.length));
		int r2 = (int)(Math.random() * ( (gene.length) - 1));
		
		while(r2 == r)
			r2++;
		
		int t = gene[r];
		gene[r] = gene[r2];
		gene[r2] = t;
	}

	private static void crossover(int[] gene1, int[] gene2)
	{
		int r = ( (int)(Math.random() * ((gene1.length) - 2)) ) + 1;
		
		int[] t = new int[gene1.length]; int[] t2 = new int[gene1.length];
		
		System.arraycopy(gene1, 0, t2, 0, r);
		System.arraycopy(gene2, r, t2, r, ((gene1.length) - r));
		
		System.arraycopy(gene2, 0, t, 0, r);
		System.arraycopy(gene1, r, t, r, ((gene1.length) - r));
		
		System.arraycopy(t, 0, gene1, 0, gene1.length);
		System.arraycopy(t2, 0, gene2, 0, gene1.length);
		
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

	private static boolean doCross(double p1_x, double p1_y, double p2_x, double p2_y, double p3_x, double p3_y, double p4_x, double p4_y)
	{
		int a = orientation(p1_x, p1_y, p2_x, p2_y, p3_x, p3_y);
		int b = orientation(p1_x, p1_y, p2_x, p2_y, p4_x, p4_y);
		int c = orientation(p3_x, p3_y, p4_x, p4_y, p1_x, p1_y);
		int d = orientation(p3_x, p3_y, p4_x, p4_y, p2_x, p2_y);

		return ( (a != b) && (c != d) );
	}

	//p, q, r
	private static int orientation(double p1_x, double p1_y, double p2_x, double p2_y, double p3_x, double p3_y) 
	{
	    double val = (p2_y - p1_y) * (p3_x - p2_x) - (p2_x - p1_x) * (p3_y - p2_y); 
	  	//System.out.println("************val = " + val);System.out.printf("************%f|%f - %f|%f - %f|%f\n", p1_x, p1_y, p2_x, p2_y, p3_x, p3_y);
	    if (val == 0.0) return 0; // colinear 
	  
	    return (val > 0.0) ? 1 : 2; // clockwise / counterclockwise 
	} 
	
	private static int countCrossings(int ordering[])
	{
		int ret = 0; int mov = 200;
		int radius = 100; double chunk = ( (Math.PI * 2) / ( (double)N ) );
		int h = 0;
		for(int i = 0; i < ordering.length; i++)
		{
			for(int j = i + 1; j < ordering.length; j++)
			{
				if( adjacencyMatrix[ ordering[i] ][ ordering[j] ] == 1 )
				{
					for(int k = i + 1; k < ordering.length; k++)
					{
						for(int m = k + 1; m < ordering.length; m++)
						{
							if( adjacencyMatrix[ ordering[k] ][ ordering[m] ] == 1 && k != j && m != j )
							{
								//System.out.println("Check lines {" + ordering[i] + ", " + ordering[j] + "} and {" + ordering[k] + ", "
								//	+ ordering[m] + "}\n");

								double vi_x = ( Math.cos(i * chunk) * radius ) + mov;
								double vi_y = ( Math.sin(i * chunk) * radius ) + mov;
								double vj_x = ( Math.cos(j * chunk) * radius ) + mov;
								double vj_y = ( Math.sin(j * chunk) * radius ) + mov;
								double vk_x = ( Math.cos(k * chunk) * radius ) + mov;
								double vk_y = ( Math.sin(k * chunk) * radius ) + mov;
								double vm_x = ( Math.cos(m * chunk) * radius ) + mov;
								double vm_y = ( Math.sin(m * chunk) * radius ) + mov;
								
								if( doCross(vi_x, vi_y, vj_x, vj_y, vk_x, vk_y, vm_x, vm_y) )
									ret++;
							}
						}	
					}
				}
			}
		}

		return ret;
	}

	private static void getBetweenessCentralities()
	{
		int closed_and_open_list[] = new int[N];

		ArrayList< ArrayList<Integer> > paths = new ArrayList< ArrayList<Integer> >();

		for(int i = 0; i < N; i++)
		{
			for(int j = i + 1; j < N; j++)
			{
				if(adjacencyMatrix[i][j] == 1)
					continue;

				for(int k = 0; k < N; k++)
				{
					if(adjacencyMatrix[i][k] == 1)
					{
						ArrayList<Integer> temp = new  ArrayList<Integer>();
						temp.add(i); temp.add(k);
						paths.add(temp);
						closed_and_open_list[k] = 1;
					}
				}

				closed_and_open_list[i] = 1;

				populatePaths(j, closed_and_open_list, paths);

				paths.clear();

				Arrays.fill(closed_and_open_list, 0);
			}
		}
	}

	@SuppressWarnings("deprecation")
	private static void populatePaths(int goal, int[] closed_and_open_list, ArrayList< ArrayList<Integer> > paths)
	{
		int times_passed_through[] = new int[N];

		boolean goal_reached = false;

		while( ( contains(closed_and_open_list, 0) ) && ( !goal_reached ) )
		{
			int size = paths.size();

			break_here:
			for(int i = 0; i < size; i++)
			{
				int path_length = paths.get(i).size();

				if(adjacencyMatrix[goal][paths.get(i).get(paths.get(i).size() - 1)] == 1)
				{
					paths.get(i).add(goal);
					goal_reached = true;
					continue break_here;
				}

				boolean branch = false;
				for(int j = 0; j < N; j++)
				{
					if(adjacencyMatrix[j][paths.get(i).get(path_length - 1)] == 1
						&& closed_and_open_list[j] == 0)
					{
						if(branch)
						{
							ArrayList<Integer> temp = new ArrayList<Integer>();

							for(int k = 0; k < path_length; k++)
							{
								temp.add(new Integer(paths.get(i).get(k)));
							}

							temp.add(j);

							paths.add(temp);
						}
						else
						{
							paths.get(i).add(j);
							branch = true;
						}
					}
				}
			}

			for(int i = 0; i < paths.size(); i++)
			{
				closed_and_open_list[ paths.get(i).get(paths.get(i).size() - 1) ] = 1;

				if( paths.get(i).get(paths.get(i).size() - 1) == goal)
					goal_reached = true;
			}
		}
		
		int count = 0;
		while(count < paths.size())
		{
			if( paths.get(count).get(paths.get(count).size() - 1) == goal)
				count++;
			else
				paths.remove(count);
		}
		
		int min = N;
		for(int i = 0; i < paths.size(); i++)
			if(paths.get(i).size() < min)
				min = paths.get(i).size();
		
		count = 0;
		while(count < paths.size())
		{
			if( paths.get(count).size() > min)
				paths.remove(count);
			else
				count++;
		}
		
		for(int i = 0; i < paths.size(); i++)
			for(int j = 1; j < paths.get(i).size() - 1; j++)
				times_passed_through[ paths.get(i).get(j) ]++;

		for(int i = 0; i < times_passed_through.length; i++)
			betweenessCentralities[i] += (double)( (double)times_passed_through[i] ) / ( (double)paths.size() );
	}
	
	private static boolean contains(int[] arr, int check)
	{
		for(int i : arr)
			if(i == check)
				return true;

		return false;	
	}

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