from abc import ABC, abstractmethod
import os
import json
from collections import Counter
from csv import writer
import pandas as pd
import matplotlib.pyplot as plt
import click


'''
Python script to generate the histogram plot from data.fasta file
'''
class DataLoader(ABC):

    @abstractmethod
    def load_data(self, filepath, file_name):
        pass

class CSVLoader(DataLoader):

    # method to load data.fasta file as input and generate the graph
    def load_data(filepath, file_name):

        # file_name = input("enter extracted file_name.txt: ")
        if os.path.isfile(file_name):
            click.echo(f"{file_name} file already exists:)")
            pass
        else:
            click.echo(f"writing your {file_name}, wait...")
            with open(filepath, 'r') as file_obj:
                file_data = file_obj.readlines()

                sequence_dict = {}

                for sequence in file_data:
                    if sequence[0] == '>':
                        continue 
                    else:
                        if sequence not in sequence_dict:
                            sequence_dict.update({f'{sequence.rstrip()}':file_data.count(sequence)})
                        else:
                            pass

            with open(file_name, "w") as sequence_file:
                sequence_file.write(json.dumps(sequence_dict))
            click.echo("Your txt file is ready to view:)")
        with open(file_name) as sequence_obj:
            textline = sequence_obj.readlines()[0]

        return json.loads(textline)
    
    def write_csv(obj_dict, csv_file):

        if os.path.isfile(csv_file):
            # click.echo("your lengthwise occurrence csv file exists:)")
            click.echo(f"your {csv_file} csv file exists:)")
            pass
        else:
            type_csv = click.prompt("type of your csv: lengthwise_occurrence/sequence_occurrence")
            h1 = click.prompt(f"give {csv_file} column1 header")
            h2 = click.prompt(f"give {csv_file} column2 header")

            def lengthwise_occurrence():
                # csv_file = input("Enter your lengthwise occurrence csv file in format any file_name.csv: ")
 
                # creating lengthwise unique sequence dictionary
                c = Counter()
                for key,value in obj_dict.items():
                    if key in obj_dict.keys():
                        c[len(key)]+=value

                c.items()

                # creating two lists for two columns
                sequence_lengthwise = []
                lengthwise_count = []
                for key,values in c.items():
                    sequence_lengthwise.append(key)
                    lengthwise_count.append(values)
                return sequence_lengthwise, lengthwise_count

            def sequence_occurrence():
                sequence_names = list(obj_dict.keys())
                sequence_count = list(obj_dict.values())
                return sequence_names, sequence_count

            click.echo(f"creating your {csv_file} csv file...")
            if type_csv == "sequence_occurrence":
                col1, col2 = sequence_occurrence()
            else:
                col1, col2 = lengthwise_occurrence()

            with open(csv_file, "w", encoding="utf-8", newline="") as file:
                thewriter = writer(file)
                header = [h1, h2]
                thewriter.writerow(header)
                for item in zip(col1, col2):
                    thewriter.writerow(item)

            # click.echo(f"creating percentages of sequence_occurrences and writing in {csv_file}")
            df = pd.read_csv(csv_file)
            count = Counter(dict(zip(df[h1],df[h2])))
            # count
            total_sequence_distribution = count.total()
            # total_sequence_distribution
            # sequences = []
            percentages = []
            for key,val in count.items():
                # sequences.append(key)
                percentages.append((val*100/total_sequence_distribution))
            
            h3 = click.prompt(f"give {csv_file} column3 header for percentages of sequence_occurrences")
            df[h3] = percentages
            df.to_csv(csv_file, index=False)
            # click.echo(f"added percentages version of the sequence_occurrence column to {csv_file}.")
            click.echo(f"Your {csv_file} is ready to view:)")


    def lengthwise_histogram(csv_file):
        # read the csv file
        # df = pd.read_csv('lengthwise_occur1.csv')
        df = pd.read_csv(csv_file)
        # print(df)

        # create histogram
        # figsize = click.prompt("Enter figure size numbers as tuple", type=int, default=(10,5))
        # l = click.prompt("Enter figure size l,w as tuple", type=click.Choice([str(i) for i in range(1,21)]), show_choices=True, default=10) #choices are strings, so input have to string
        # w = click.prompt("Enter figure size l,w as tuple", type=click.Choice([str(i) for i in range(1,21)]), show_choices=True, default=5)
        # l = int(l)
        # w = int(w)
        l = int(click.prompt("Enter figure size l,w as tuple", type=click.Choice([str(i) for i in range(1,21)]), show_choices=True, default=str(10))) #choices are strings, so input have to string
        w = int(click.prompt("Enter figure size l,w as tuple", type=click.Choice([str(i) for i in range(1,21)]), show_choices=True, default=str(5)))
        # print(figsize)
        # for l,w in figsize:
        # l = int(l)
        # w = int(w)

        # figure, axes = plt.subplots(1, 1, figsize=(10,5))
        figure, axes = plt.subplots(1, 1, figsize=(l,w))
        # figure, axes = plt.subplots
        x_range_min = int(click.prompt("x axis minimum range", default=df['sequence_lengthwise'].min()))
        x_range_max = int(click.prompt("x axis maximum range", default=df['sequence_lengthwise'].max()))
        # y_range_min = click.prompt("y axis minimum range", default=1)
        y_range_max = int(click.prompt("y axis maximum range", default=df['lengthwise_count'].max()))
        bins = int(click.prompt("bins for the hist2d", default=df['sequence_lengthwise'].max()))

        histogram = axes.hist2d(df['sequence_lengthwise'],df['lengthwise_count'], range=[[x_range_min, x_range_max], [1,y_range_max]], bins=bins)
        # histogram = axes.hist2d(df['sequence_lengthwise'],df['lengthwise_count'], range=[[15, 65], [1,200000]], bins=65)#, norm=LogNorm())
        # histogram = axes.hist2d(df['sequence_lengthwise'],df['lengthwise_count'], range=[[df['sequence_lengthwise'].min(), df['sequence_lengthwise'].max()], [1,df['lengthwise_count'].max()]], bins=df['sequence_lengthwise'].max())
        # histogram = axes.hist2d(df['sequence_lengthwise'],df['lengthwise_count'], range=[[df['sequence_lengthwise'].min(), df['sequence_lengthwise'].max()], [1,df['lengthwise_count'].max()]], bins=65)
        # histogram = axes.hist2d(df['sequence_lengthwise'],df['lengthwise_count'], range=[[df['sequence_lengthwise'].min(), df['sequence_lengthwise'].max()], [df['lengthwise_count'].min(),df['lengthwise_count'].max()]], bins=65)

        x_label = click.prompt("Label for x axis", default="sequences_lengthwise")
        # plt.xlabel('sequences_lengthwise')
        plt.xlabel(x_label)
        y_label = click.prompt("Label for y axis", default="lengthwise_count")
        # plt.ylabel('lengthwise_count')
        plt.ylabel(y_label)
        plot_title = click.prompt("Label for graph title", default="Lengthwise frequency distribution")
        # plt.title("Lengthwise frequency distribution")
        plt.title(plot_title)

        # create colorbar
        cbar = figure.colorbar(histogram[3])
        cbar_label = click.prompt("Label for vertical color bar", default="occurrence bar")
        size = click.prompt("size of the label", default=10)
        # cbar.set_label(label="occurrence bar", size=10)
        cbar.set_label(label=cbar_label, size=size)
        
        figure_name = click.prompt("Name of your figure to be saved", default="histo.png")
        facecolor = click.prompt("facecolor of the figure", default="white")
        figure.savefig(figure_name, facecolor=facecolor)
        # figure.savefig('lengthwise.png', facecolor='white') # Note: for default python py -3.11 giving filepath as only "lengthwise.png" getting image saved in SMcodes7 instead of folder because matplotlib can access active main directory, so provide Job_task_matplotlib/lengthwise.png path to the saved image
        click.echo("Your histogram is ready to view:)")
        # showing the histogram
        plt.show()
        plt.close(figure)


    def barh_maximum_sequences(csv_file):
        # reading csv
        colbhx = click.prompt("give your csv file column name for barh x-axis")
        colbhy = click.prompt("give your csv file column name for barh y-axis")
        df = pd.read_csv(csv_file, usecols=[colbhx, colbhy])

        max = click.prompt("Enter number of maximum sequences to plot", default=10)
        maxn_sequences = df.nlargest(max, [colbhy])

        # max_sequences
        x = [f'{seq} ({len(seq)})' for seq in maxn_sequences[colbhx]]
        y = maxn_sequences[colbhy]

        l = int(click.prompt("Enter figure size l,w as tuple", type=click.Choice([str(i) for i in range(1,40)]), show_choices=True, default=str(31))) #choices are strings, so input have to string
        w = int(click.prompt("Enter figure size l,w as tuple", type=click.Choice([str(i) for i in range(1,40)]), show_choices=True, default=str(10)))
        figure, axes = plt.subplots(1,1, figsize=(l,w))#(31,10))#(38,10)) 

        # plot barh
        bar_color = click.prompt("bars color", default="purple")
        bargraph = axes.barh(x,y, color=bar_color)

        # Add labels and padding between axes and tick labels
        xtick_labelsize = click.prompt("give barh x axis ticks labelsize", default=10)
        xtick_labelcolor = click.prompt("give barh x axis tick labelcolor", default='magenta')
        xtick_padding = click.prompt("give barh x axis ticks padding size", default=5)
        axes.xaxis.set_tick_params(labelsize=xtick_labelsize, labelcolor=xtick_labelcolor, pad=xtick_padding) # labelsize=20 works # width=10 increases scale ticks to black bars
        
        ytick_labelsize = click.prompt("give barh y axis labelsize", default=10)
        ytick_labelcolor = click.prompt("give barh y axis tick labelcolor", default='magenta')
        ytick_padding = click.prompt("give barh y axis ticks padding size", default=10)
        axes.yaxis.set_tick_params(labelsize=ytick_labelsize, labelcolor=ytick_labelcolor, pad=ytick_padding) # size=50, inceases scale ticks to long lines

        # Add labels # plt y is barh x and plt x is barh y
        labelx = click.prompt("give barh x axis label", default="Unique Sequences and their lengths")
        font_sizex = click.prompt("give barh x axis fontsize", default=15)
        labelx_color = click.prompt("give barh x axis label color", default='purple')
        plt.ylabel(labelx, fontsize=font_sizex, color=labelx_color, fontweight='bold')

        labely = click.prompt("give barh y axis label", default="Sequence_occurence")
        font_sizey = click.prompt("give barh y axis fontsize", default=15)
        labely_color = click.prompt("give barh y axis label color", default='purple')
        plt.xlabel(labely, fontsize=font_sizey, color=labely_color, fontweight='bold')

        # graph title
        graph_title = click.prompt("give title to your barh graph", default="10 Sequences with maximum frequency of distribution")
        pad = click.prompt("padding from graph", default=10)
        title_size = click.prompt("give font size", default=15)
        title_color = click.prompt("title color", default='purple')
        plt.title(graph_title, pad=pad, fontsize=title_size, fontweight='bold', color=title_color)

        # Add annotation to bars
        bar_top_textcolor = click.prompt("color of text on bar tops", default='magenta')
        bar_top_textsize = click.prompt("fontsize of text on bar tops", default=12)
        for i in axes.patches:
            plt.text(i.get_width()+0.2, i.get_y()+0.5,
                    str(round((i.get_width()), 2)),
                    fontsize=bar_top_textsize, fontweight='bold',
                    color=bar_top_textcolor)
        
        # save te figure
        figure_file = click.prompt("give your figure name as .png/jpg/svg", default='barh.png')
        figure.savefig(figure_file, facecolor='white')
        plt.show()
        # close the figure
        plt.close(figure)

    def pie_chart_maximum_sequences(csv_file):
        
        # reading csv
        col_pie_names = click.prompt("give your csv file column name for sequences names or legend display list")
        col_pie_percentages = click.prompt("give your csv file column name for percentages")
        df = pd.read_csv(csv_file, usecols=[col_pie_names, col_pie_percentages])

        # max10_percentages = df['Percentages'].nlargest(10,keep='all') # returns the column with index positions
        # max10_percentages # correct 1st method
        # max10_percentages = df.nlargest(10, ['Percentages'], keep='all') # returns whole df table with index positions
        # max10_percentages # correct 2nd method
        max = click.prompt("Enter number of maximum sequences to plot", default=10)
        maxn_percentages = df.nlargest(max, [col_pie_percentages], keep="all")

        # Note: can also use zip for two lists
        # list(zip(list1,list2))
        labels = [f'{seq} ({len(seq)})' for seq in maxn_percentages[col_pie_names]]
        pies = maxn_percentages[col_pie_percentages]
        
        # creating pie figure
        l = int(click.prompt("Enter figure size l,w as tuple", type=click.Choice([str(i) for i in range(1,40)]), show_choices=True, default=str(15))) #choices are strings, so input have to string
        w = int(click.prompt("Enter figure size l,w as tuple", type=click.Choice([str(i) for i in range(1,40)]), show_choices=True, default=str(8)))
        fig, ax = plt.subplots(figsize=(l,w))
        # explode = (0.1, 0.1, 0.1, 0.1, 0, 0, 0, 0, 0, 0)
        total = sum(maxn_percentages[col_pie_percentages])
        # to show absolute values from the pies list instead of relative percentages
        wedges, texts, autotexts = ax.pie(pies, autopct=lambda p: '{:.0f}%'.format(p * total / 100), normalize=True) #,autopct='%.1f%%') #,textprops={"fontsize":10} #, labels=labels) #frame=True) # Pass a function or format string to autopct to label slices
        legend_title = click.prompt("give title to the legend", default="Top 10 Maximum Occurring Sequences with their length")
        font_size = click.prompt("give legend lists's fontsize", default=10)
        mylegend = ax.legend(wedges, labels, title=legend_title, fontsize=font_size, loc='lower center', bbox_to_anchor=(1,0.5,0.5,1)) #ncol=2

        # change plt's title
        autotext_size = click.prompt("give size to the texts inside the pies", default=10)
        plt.setp(autotexts, size=autotext_size, weight="bold") # this is making pies percentage texts bold

        chart_title = click.prompt("give title to your pie chart", default="Pie-chart: Top ten maximum occurring Sequences (in %)")
        title_size = click.prompt("give size to pie chart title", default=20)
        ax.set_title(chart_title, fontsize=title_size) #, fontweight='bold')

        # save your figure
        figure_name = click.prompt("give name to your figure file as .png/.jpg/.svg", default="pie.png")
        fig.savefig(figure_name, facecolor='white')
        plt.show()
        plt.close(fig)


@click.command(help="DNA Sequences Data Analysis tool using fasta, txt, and csv files data to create histogram, barh and pie chart plots")
@click.option("-n", "--files", multiple=True, type=str, help="Files") #type=click.File()
@click.option("-t", "--type", "command", type=click.Choice(["fasta", "csv", "histo", "barh", 'pie'], case_sensitive=False), default="fasta", help="Type of operation")
def main(files, command):
    if command == "pie":
        CSVLoader.pie_chart_maximum_sequences(files[0])
    elif command == "barh":
        CSVLoader.barh_maximum_sequences(files[0])
    elif command == "histo":
       CSVLoader.lengthwise_histogram(files[0]) # giving single csv file input works
    elif command == "csv":
        ctx = main.callback(files, command="fasta") # but needs two file values
        # print(ctx)
        # print(files) # if given single value it's a single value tuple
        CSVLoader.write_csv(ctx, files[2]) # so give three file values even if not needed
    elif command == "fasta":
        # ctx = files # is a tuple due to multiple =True
        ctx = CSVLoader.load_data(files[0],files[1])
        return ctx

if __name__ == '__main__':
    main()
