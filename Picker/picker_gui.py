import os
import math
import obspy
import numpy as np
import tkinter as tk
import matplotlib.pyplot as plt

from math import pi
from math import sin,cos,tan
from tkinter import messagebox
from tkinter import filedialog
from obspy.taup import TauPyModel
from obspy.imaging.mopad_wrapper import beach
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


"""
    P-wave Polarity Picker

    Author: Fu Yin
    Date: 2023-03-06
    Email: fy21@rice.edu

"""


def get_taup_tp_ts(model,depth,distance,degree=None):
    if degree==False:
        distance = distance/111.19

    time_p = model.get_travel_times(source_depth_in_km=depth,
                                distance_in_degree=distance,
                                phase_list=["p", "P"])

    time_s = model.get_travel_times(source_depth_in_km=depth,
                                distance_in_degree=distance,
                                phase_list=["s", "S"])

    ray_p = time_p[0].ray_param
    tp = time_p[0].time
    angle_p = time_p[0].incident_angle

    ray_s = time_s[0].ray_param
    ts = time_s[0].time
    angle_s = time_s[0].incident_angle

    return ray_p,tp,angle_p,ray_s,ts,angle_s


def project_beachball(AZM, TKO, R=1, menthod='schmidt'):
    AZM = AZM/180*pi
    TKO = TKO/180*pi

    # Schmidt (Lambert, equal-area) default
    if menthod=='schmidt':
        r = math.sqrt(2)*sin(TKO/2)
    # Wulff projection (Stereographic, equal-angle) not recommmended
    elif menthod=='wulff':
        r = tan(TKO/2)
    else:
        raise ValueError('projection error!')

    X = R*r*sin(AZM)+R
    Y = R*r*cos(AZM)+R

    return X, Y


class GraphView:
    def __init__(self, master):
        self.master = master
        self.init_plot()
        self.info = []
        self.signs = []
        self.markers = []
        self.is_pressed = False
        self.x0, self.y0 = None, None
        self.fig.canvas.mpl_connect('pick_event', self.on_pick)
        self.fig.canvas.mpl_connect('scroll_event', self.scroll_triger)
        self.fig.canvas.mpl_connect('button_press_event', self.on_button_press)
        self.fig.canvas.mpl_connect('button_release_event', self.on_button_release)
        self.fig.canvas.mpl_connect('motion_notify_event', self.on_motion)

        
    def init_plot(self):
        self.fig, self.ax = plt.subplots(1,2,gridspec_kw=dict(width_ratios=[2, 1]),figsize=(12, 4), dpi=100)
        beach1 = beach(np.array([0,90,0]), xy=(50, 50),  linewidth=1,width=100-1, alpha=1,\
        facecolor='g',bgcolor='w', edgecolor='k',mopad_basis='NED',nofill=False,zorder=1 )
        self.ax[1].add_collection(beach1) 
        self.ax[1].set_aspect("equal")
        self.ax[1].set_axis_off() 
        self.ax[1].set_xlim(0,100)
        self.ax[1].set_ylim(0,100)
        self.ax[1].set_title("Strike=%.1f, Dip=%.1f, Rake=%.1f" % (0, 90, 0))
        self.ax[0].set_xlabel('Time(s)')
        self.ax[0].set_ylabel('Trace')
        self.fig.suptitle('P-wave Polarity Picker')


    def XY_init(self):
        self.x_start = 0
        self.x_end = self.max_time
        self.y_start = -1
        self.y_end = self.ydist
        self.ax[0].set_xlim(self.x_start, self.x_end)
        self.ax[0].set_ylim(self.y_start, self.y_end)
        self.fig.canvas.draw()


    def clear(self):
        self.ax[0].clear()
        self.ax[1].clear()
        beach1 = beach(np.array([0,90,0]), xy=(50, 50),  linewidth=1,width=100-1, alpha=1,\
        facecolor='g',bgcolor='w', edgecolor='k',mopad_basis='NED',nofill=False,zorder=1 )
        self.ax[1].add_collection(beach1) 
        self.ax[1].set_aspect("equal")
        self.ax[1].set_axis_off() 
        self.ax[1].set_xlim(0,100)
        self.ax[1].set_ylim(0,100)
        self.ax[1].set_title("Strike=%.1f, Dip=%.1f, Rake=%.1f" % (0, 90, 0))
        self.ax[0].set_xlabel('Time(s)')
        self.ax[0].set_ylabel('Trace')
        self.fig.suptitle('P-wave Polarity Picker')
        self.info = []
        self.markers = []
        self.signs = []
        self.is_pressed = False
        self.x0, self.y0 = None, None
        self.fig.canvas.draw()


    def on_button_press(self, event):
        if event.inaxes != self.ax[0]:
            return
        if event.button == 1:
            self.is_pressed = True
            self.x0, self.y0 = event.xdata, event.ydata


    def on_button_release(self, event):
        if event.button == 1:
            self.is_pressed = False


    def on_motion(self, event):
        if not self.is_pressed:
            return
        dx = event.xdata - self.x0
        dy = event.ydata - self.y0
        xlim = self.ax[0].get_xlim()
        ylim = self.ax[0].get_ylim()
        self.x_start = xlim[0] - dx
        self.x_end = xlim[1] - dx
        self.y_start = ylim[0] - dy
        self.y_end = ylim[1] - dy
        self.ax[0].set_xlim(xlim - dx)
        self.ax[0].set_ylim(ylim - dy)
        self.fig.canvas.draw()


    def scroll_triger(self, event):
        x_mid, y_mid = event.xdata, event.ydata
        x_len_r = self.x_end - x_mid
        x_len_l = x_mid - self.x_start
        y_len_r = self.y_end - y_mid
        y_len_l = y_mid - self.y_start
        scale = 0.8
        if event.button == 'up':
            self.x_start = x_mid - x_len_l*scale 
            self.x_end = x_mid + x_len_r*scale
            self.y_start = y_mid - y_len_l*scale
            self.y_end = y_mid + y_len_r*scale
            self.ax[0].set_xlim(self.x_start, self.x_end)
            self.ax[0].set_ylim(self.y_start, self.y_end)
        if event.button == 'down':
            self.x_start = x_mid - x_len_l/scale
            self.x_end = x_mid + x_len_r/scale
            self.y_start = y_mid - y_len_l/scale
            self.y_end = y_mid + y_len_r/scale
            self.ax[0].set_xlim(self.x_start, self.x_end)
            self.ax[0].set_ylim(self.y_start, self.y_end)
        for i in range(0,len(self.ax[0].texts)):
            self.ax[0].texts[i].set_clip_on(True)
        self.fig.canvas.draw()


    def plot_waveform(self, tr):
        self.ydist = 0
        self.max_time = 0
        for i in range(0, len(tr)):
            dt = tr[i].stats.delta
            npts =tr[i].stats.npts
            time = np.arange(0, npts*dt, dt)
            if max(time) > self.max_time:
                self.max_time = max(time)
            data = tr[i].data/np.max(tr[i].data)/2
            name = tr[i].stats.network + '_' + tr[i].stats.station + '_' + tr[i].stats.channel
            self.info.append(name)
            self.ax[0].plot(time, data+self.ydist, 'b-', picker=5, alpha = 0.8)
            self.ax[0].text(1, self.ydist+0.25, name, va='center', ha='left')
            self.ydist += 1
        self.x_start = 0
        self.x_end = self.max_time
        self.y_start = -1
        self.y_end = self.ydist
        self.ax[0].set_xlim(0,self.max_time)
        self.ax[0].set_ylim(-1,self.ydist)


    def plot_beachball(self, FM, X, Y, Name):
        self.X = X
        self.Y = Y
        self.Name = Name
        beach1 = beach(FM, xy=(50, 50),  linewidth=1,width=100-1, alpha=1,\
                facecolor='g',bgcolor='w', edgecolor='k',mopad_basis='NED',nofill=False,zorder=1 )
        self.ax[1].add_collection(beach1) 
        self.ax[1].set_aspect("equal")

        for i in range(0, len(X)):
            self.ax[1].plot(X[i], Y[i], "rv", ms=10,zorder=1) 
            self.ax[1].text(X[i], Y[i], Name[i], horizontalalignment='right', verticalalignment='center',\
                fontsize=10, color='black',bbox = dict(facecolor = "r", alpha = 0.0),zorder=1) 

        self.ax[1].set_xlim(0,100)
        self.ax[1].set_ylim(0,100)
        self.ax[1].set_title("Strike=%.1f, Dip=%.1f, Rake=%.1f" % (FM[0], FM[1], FM[2]))



    def on_pick(self, event):
        x, y = event.artist.get_data()
        ind = event.ind[0]
        t, a = x[ind], y[ind]
        
        if event.mouseevent.button == 1:
            color = 'r'  
            ss = '+'
        elif event.mouseevent.button == 3:
            color = 'k'
            ss = '-'
        else:
            return

        marker = self.ax[0].plot(t, a, color+'o', markersize=5)[0]
        i = round(a)
        sign = self.ax[1].text(self.X[i], self.Y[i], ss, horizontalalignment='left', verticalalignment='center',\
                fontsize=15, color='black',bbox = dict(facecolor = color, alpha = 0.3),zorder=2) 
        
        self.signs.append(sign)
        self.markers.append(marker)
        self.fig.canvas.draw()


    def delete_last_pick(self):
        if self.markers:
            marker = self.markers.pop()
            marker.remove()
            self.fig.canvas.draw()

        if self.signs:
            sign = self.signs.pop()
            sign.remove()
            self.fig.canvas.draw()


    def save_picks(self, filename, AZM_all, TKO_all):
        with open(filename, 'w') as f:
            f.write("Trace_Num\tArrive_Time\tNetwork_Station_Channel\tPolarity\tAzimuth\tTakeoff\n")
            for marker in self.markers:
                pos = marker.get_data()
                x = pos[0][0]
                y = pos[1][0]
                color = marker.get_color()[0]
                if color == 'r':
                    polarity = '+'
                elif color == 'k':
                    polarity = '-'
                else:
                    raise ValueError("Unknown color")
                name = self.info[int(y)]
                azm = AZM_all[int(y)]
                tko = TKO_all[int(y)]
                f.write(f"{round(y):d}\t{x:f}\t{name:s}\t{polarity}\t{azm:.2f}\t{tko:.2f}\n")
        messagebox.showinfo("Save Successfully", f"The data is saved into '{filename}'")



class MainWindow:
    def __init__(self, master):
        self.master = master
        self.master.geometry("800x600")
        self.master.title("Picker")
        self.graph_view = GraphView(self.master)
        
        frame_4 = tk.Frame(self.master)
        frame_4.pack(side=tk.BOTTOM)
        self.vmodel_label = tk.Label(frame_4, text='V_model')
        self.vmodel_label.pack(side="left")
        self.vmodel_entry = tk.Entry(frame_4, width=14)
        self.vmodel_entry.pack(side="left")

        frame_3 = tk.Frame(self.master)
        frame_3.pack(side=tk.BOTTOM)
        self.rake_label = tk.Label(frame_3, text=' Rake')
        self.rake_label.pack(side="left")
        self.rake_entry = tk.Entry(frame_3, width=8)
        self.rake_entry.pack(side="left")

        frame_2 = tk.Frame(self.master)
        frame_2.pack(side=tk.BOTTOM)
        self.dip_label = tk.Label(frame_2, text='   Dip')
        self.dip_label.pack(side="left")
        self.dip_entry = tk.Entry(frame_2, width=8)
        self.dip_entry.pack(side="left")

        frame_1 = tk.Frame(self.master)
        frame_1.pack(side=tk.BOTTOM)
        self.strike_label = tk.Label(frame_1, text='Strike')
        self.strike_label.pack(side="left")
        self.strike_entry = tk.Entry(frame_1, width=8)
        self.strike_entry.pack(side="left")

        button_frame = tk.Frame(master)
        button_frame.pack(side=tk.BOTTOM, padx=10, pady=10)
        self.load_button = tk.Button(button_frame, text="Load", command=self.load_data)
        self.load_button.pack(side="left")
        self.delete_button = tk.Button(button_frame, text="Delete", command=self.on_delete)
        self.delete_button.pack(side="left")
        self.save_button = tk.Button(button_frame, text="Save", command=self.on_save)
        self.save_button.pack(side="left")
        self.plot_button = tk.Button(button_frame, text="Plot", command=self.plot)
        self.plot_button.pack(side="left")
        self.init_button = tk.Button(button_frame, text="XY_init", command=self.XY_init)
        self.init_button.pack(side="left")


    # recover to the initial figure size
    def XY_init(self):
        self.graph_view.XY_init()


    # load data
    def load_data(self):
        self.graph_view.clear()
        self.strike_entry.delete(0, tk.END)
        self.dip_entry.delete(0, tk.END)
        self.rake_entry.delete(0, tk.END)

        file_paths = filedialog.askopenfilenames(initialdir='./')
        if not file_paths:
            return
        
        self.tr = obspy.read(file_paths[0])
        if len(file_paths) > 1:
            for i in range(0, len(file_paths)):
                self.tr += obspy.read(file_paths[i])

        self.graph_view.plot_waveform(self.tr)
        self.graph_view.fig.canvas.draw()


    # save data
    def on_save(self):
        filename = "picks.txt"
        self.graph_view.save_picks(filename, self.AZM_all, self.TKO_all)


    # delete data
    def on_delete(self):
        self.graph_view.delete_last_pick()


    # plot data
    def plot(self):
        v_model = self.vmodel_entry.get()
        model_path = os.path.join(v_model) 
        if v_model == 'iasp91':
            model = TauPyModel(model='iasp91') 
        elif v_model == 'prem':
            model = TauPyModel(model='prem')
        else:
            model = TauPyModel(model=model_path)

        for i in range(0,len(self.tr),1):
            depth = self.tr[i].stats.sac['evdp']
            distance = self.tr[i].stats.sac['dist']
            ray_p,tp,angle_p,ray_s,ts,angle_s = get_taup_tp_ts(model,depth,distance,degree=False)                                            
            self.tr[i].stats.sac["user1"]=angle_p
            self.tr[i].stats.sac["user2"]=angle_s

        menthod='schmidt'   # 'schmidt' 'wulff'
        X=[]; Y=[]; Name=[]; 
        self.AZM_all=[]; self.TKO_all=[]
        for i in range(0,len(self.tr)):
            AZM = self.tr[i].stats.sac['az'] 
            self.AZM_all.append(AZM)
            TKO = self.tr[i].stats.sac['user1']
            self.TKO_all.append(TKO)
            name = self.tr[i].stats.network+'_'+self.tr[i].stats.station    
            x, y = project_beachball(AZM, TKO, R=100/2, menthod=menthod)   
            X.append(x) 
            Y.append(y)
            Name.append(name)

        dip = self.dip_entry.get()
        strike = self.strike_entry.get()
        rake = self.rake_entry.get()
        if strike == '' or dip == '' or rake == '':
            FM = [0, 90, 0]
        else:
            FM = [float(strike), float(dip), float(rake)]
        self.graph_view.plot_beachball(FM, np.array(X), np.array(Y), np.array(Name))
        self.graph_view.fig.canvas.draw()



if __name__ == '__main__':
    root = tk.Tk()
    app = MainWindow(root)
    canvas = FigureCanvasTkAgg(app.graph_view.fig, master=root)
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
    canvas.draw()
    root.mainloop()