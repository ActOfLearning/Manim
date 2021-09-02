from sys import base_exec_prefix
from manim import *
from manim.opengl import *
import os
from manim.mobject.svg.style_utils import SVG_DEFAULT_ATTRIBUTES

from manim.utils.file_ops import add_version_before_extension
from numpy import arange, diff, dot
from scipy.integrate import odeint

#config.use_opengl_renderer = True

class Anglet(VGroup):

    CONFIG = {
        'radius': 0.5,
        'color': RED,
        'opacity': 0.4,
        'stroke_width': 5,
    }

    def __init__(self, A, O, B, radius = 0.5, color = RED, opacity = 0.4, stroke_width = 5, below_180 = True, **kwargs):

        VMobject.__init__(self, **kwargs)
        OA, OB = A-O, B-O
        if below_180:
            theta = np.angle(complex(*OA[:2])/complex(*OB[:2])) # angle of OB to OA
        else:
            theta = TAU + np.angle(complex(*OA[:2])/complex(*OB[:2]))

        if abs(abs(theta) - PI / 2) <= 1e-2:
            self.add(Elbow(angle = np.angle(complex(*OB[:2])), color = color, width = radius).shift(O))
        else:
            self.add(Sector(inner_radius = 0, outer_radius = radius, angle = theta, arc_center = O, color = color, fill_opacity = opacity, start_angle = Line(O, B).get_angle()))
            self.add(Arc(start_angle=Line(O, B).get_angle(), angle=theta, radius=radius,
                     stroke_width=stroke_width, color=color, arc_center=O))


class CycloidGeneration(Scene):
    def construct(self):
        t = ValueTracker(0)
        tlab = DecimalNumber(0)
        tlab.set_value(t.get_value(), num_decimal_places = 2)
        tlab.add_updater(lambda m: m.set_value(t.get_value()))
        tgp = VGroup(MathTex("t="),tlab)
        tgp.arrange(RIGHT)

        path = VMobject()
        dot, cir = Dot(), Circle().rotate_in_place(90 * DEGREES)
        cir.save_state()
        dot.move_to(cir.get_start())
        path.set_points_as_corners([dot.get_center(), dot.get_center()])

        def circle_update(obj):
            cir.restore()
            cir.rotate_in_place(t.get_value())
            cir.shift(t.get_value() * RIGHT)
            obj.become(cir)
            #dot.move_to(cir.get_start())
        def update_path(path):
            previous_path = path.copy()
            previous_path.add_points_as_corners([dot.get_center()])
            path.become(previous_path)
        
        path.add_updater(update_path)
        cir.add_updater(circle_update)
        dot.add_updater(lambda m: m.move_to(cir.get_start()))

        tgp.move_to(2 * UP)
        self.play(Write(cir), Write(dot), Write(tgp))
        self.add(path, cir, dot)
        self.play(
            t.animate.set_value(PI),
            rate_func = linear,
            run_time = 3
        )

        tlab.clear_updaters()
        dot.clear_updaters()
        path.remove_updater(update_path)
        cir.remove_updater(circle_update)
        self.wait(5)


class Brachis(MovingCameraScene):
    def construct(self):
        t = ValueTracker(0)
        tlab = DecimalNumber(0)
        tlab.set_value(t.get_value(), num_decimal_places = 2)
        tlab.add_updater(lambda m: m.set_value(t.get_value()))
        tgp = VGroup(MathTex("t="),tlab)
        tgp.arrange(RIGHT)

        self.camera.frame.shift(3 * DR + 1.5 * RIGHT)
        
        cycgrph = ParametricFunction(
            lambda u: np.array([u - np.sin(u), -1 + np.cos(u), 0]),
            t_range = np.array([0, 180 * DEGREES]), fill_opacity = 0).set_color(RED)
        cycgrph.scale_about_point(3, ORIGIN)
        sd, ed = Dot(), Dot()
        sd.move_to(cycgrph.get_start())
        ed.move_to(cycgrph.get_end())
        dotgrp = VGroup(sd, ed)
        self.play(Create(cycgrph), Write(dotgrp))

        pdot = Dot().move_to(sd.get_center()).set_color(GREEN)
        pdot.state_vector = np.array([60 * DEGREES, 0.0])
        self.play(Write(pdot))

        g, eps = 9.81, 1 / 1024 / 1024

        def pdotupdate(obj, dt):
            for x in range(2):
                theta, z = obj.state_vector
                thetadot = z
                print(theta, theta % (2 * PI), z)
                if abs(theta) % (2 * PI) < 1 / 1024 / 1024 or 2 * PI - abs(theta) % (2 * PI) < 1 / 1024 / 1024:
                    zdot = 64
                elif PI - abs(theta) % PI < 1 / 1024 / 1024:
                    zdot = 0
                else:
                    zdot = (g / 2 / 3 - z * z / 2) / (np.tan(theta / 2))
                obj.state_vector += dt * np.array([thetadot, zdot]) / 2
            pdot.move_to(3 * np.array([theta - np.sin(theta), -1 + np.cos(theta), 0]))
        pdot.add_updater(pdotupdate)
        self.add(pdot)

        self.wait(2)

        pdot.remove_updater(pdotupdate)
        self.wait(5)


class BrachisODE(MovingCameraScene):
    def construct(self):
        #https://sam-dolan.staff.shef.ac.uk/mas212/notebooks/ODE_Example.html
        #https://pastebin.com/2s9At2XE
        t = ValueTracker(0)
        tlab = DecimalNumber(0)
        tlab.set_value(t.get_value(), num_decimal_places = 2)
        tlab.add_updater(lambda m: m.set_value(t.get_value()))
        tgp = VGroup(MathTex("t="),tlab)
        tgp.arrange(RIGHT)

        self.camera.frame.shift(3 * DR + 1.5 * RIGHT)
        
        cycgrph = ParametricFunction(
            lambda u: np.array([u - np.sin(u), -1 + np.cos(u), 0]),
            t_range = np.array([0, 180 * DEGREES]), fill_opacity = 0).set_color(RED)
        cycgrph.scale_about_point(3, ORIGIN)
        sd, ed = Dot(), Dot()
        sd.move_to(cycgrph.get_start())
        ed.move_to(cycgrph.get_end())
        dotgrp = VGroup(sd, ed)
        self.play(Create(cycgrph), Write(dotgrp))

        pdot = Dot().move_to(sd.get_center()).set_color(GREEN)

        semi_dist = np.sqrt(6 * 6 + 9 * PI * PI) / 2
        lencircle = 0.75268050611004932993 * 2 * semi_dist
        thetp = np.arcsin(semi_dist / lencircle)
        alp = np.arctan(6 / 3 / PI)
        cen_height = lencircle * np.cos(thetp)
        circen = Dot(semi_dist * RIGHT + cen_height * UP)
        circen.rotate_about_origin(-alp)
        pcdot = Dot().move_to(sd.get_center()).set_color(RED)

        passingdots = VGroup(pdot, pcdot)
        self.play(Write(passingdots))

        inits, g = [1 * DEGREES / 512, 0.0], 9.81
        timevalues = np.linspace(0, 3 + 1 / 4, 1000)

        def diffeqn(vec, t):
            return [vec[1], (g / 2 / 3 - vec[1] * vec[1] / 2) / np.tan(vec[0] / 2)]
            #return [vec[1], -2 * vec[1] - 2 * vec[0] + np.cos(2 * t)]
        res = odeint(diffeqn, inits, timevalues)
        thetavalues = res[:, 0]

        cir_inits = [alp + thetp, 0]
        #cir_inits = 0
        def cir_diffeqn(u, t):
            return [u[1], - (g / lencircle) * np.sin(u[0])]
        cir_res = odeint(cir_diffeqn, cir_inits, timevalues)
        cir_thetavalues = cir_res[:, 0]

        cirgrph = ParametricFunction(
            lambda u: circen.get_center() + lencircle * np.cos(u) * DOWN + lencircle * np.sin(u) * LEFT, 
            t_range = np.array([alp - thetp, cir_inits[0]]), fill_opacity = 0
        ).set_color(GREEN)
        self.play(Write(cirgrph))
        
        def pdotupdater(obj):
            curr_time = t.get_value() / 4
            curr_theta = np.interp(curr_time, timevalues, thetavalues)
            curr_cir_theta = np.interp(curr_time, timevalues, cir_thetavalues)
            obj[0].move_to(3 * np.array([curr_theta - np.sin(curr_theta), -1 + np.cos(curr_theta), 0]))
            obj[1].move_to(circen.get_center() + lencircle * np.cos(curr_cir_theta) * DOWN + lencircle * np.sin(curr_cir_theta) * LEFT)
        
        passingdots.add_updater(pdotupdater)
        self.add(passingdots)

        self.play(
            t.animate.set_value(8),
            run_time = 8,
            rate_func = linear
        )
        
        self.wait()


class SphericalTautochrone(Scene):
    def construct(self):
        #https://sam-dolan.staff.shef.ac.uk/mas212/notebooks/ODE_Example.html
        #https://pastebin.com/2s9At2XE
        t = ValueTracker(0)
        tlab = DecimalNumber(0)
        tlab.set_value(t.get_value(), num_decimal_places = 2)
        tlab.add_updater(lambda m: m.set_value(t.get_value()))
        tgp = VGroup(MathTex("t="),tlab)
        tgp.arrange(RIGHT)

        total_time = 5
        g = 9.81
        lmda = 2 * g * (total_time / PI) * (total_time / PI)
        lmda = 1 / lmda
        lmda = 1 - lmda
        lmda = np.arccos(lmda)

        def phi(theta):
            ct, cl = np.cos(theta), np.cos(lmda)
            return (1 / np.tan(lmda / 2)) * np.arctanh(np.sqrt((cl - ct) / (1 + cl))) - np.arctan(np.sqrt((cl - ct) / (1 - cl)))
        
        def spherical_diff(theta, ti, start_angle):
            sa = max(start_angle, lmda)
            ct, cs = np.cos(theta), np.cos(sa - 1 / 1024 / 1024)
            return (PI / total_time) * np.sqrt((cs - ct) / (1 - ct))
        
        timevalues = np.linspace(0, total_time, 500)

        self.camera.scale(1 / 3)
        self.camera.light_source.move_to(3*IN)

        sph = Sphere(fill_opacity = 1 / 2)
        ax = ThreeDAxes()
        base_plane = NumberPlane().shift(IN)
        self.add(sph, base_plane)
        self.play(self.camera.animate.set_euler_angles(phi = 60 * DEGREES, theta = 30 * DEGREES))
        #self.play(self.camera.animate.set_euler_angles(phi = 90 * DEGREES, theta = 0 * DEGREES))
        #self.play(self.camera.animate.set_euler_angles(phi = 0 * DEGREES, theta = 0 * DEGREES))
        #self.play(self.camera.animate.set_euler_angles(phi = 90 * DEGREES, theta = 90 * DEGREES))

        sph_tauto = ParametricFunction(
            lambda u: np.array([np.sin(u) * np.cos(phi(u)), np.sin(u) * np.sin(phi(u)), np.cos(u)]), 
            t_range = np.array([lmda, PI - 1 / 1024]), fill_opacity = 0, color = RED
        )
        self.play(Write(sph_tauto))

        n = 3
        pdots = VGroup(*[Dot3D(radius = 1 / 32) for k in range(n)])
        gps = (PI - lmda) / n
        resvalues = []
        for j in range(n):
            u = lmda + j * gps
            vals = odeint(spherical_diff, u, timevalues, args = (u, ))
            tvals = vals[:, 0]
            resvalues += [tvals]
            pdots[j].move_to(np.sin(u) * np.cos(phi(u)) * RIGHT + np.sin(u) * np.sin(phi(u)) * UP + np.cos(u) * OUT)
        
        self.play(
            Write(pdots),
        )

        def pdotsupdater(obj):
            curr_time = t.get_value() / 4
            for j in range(n):
                u = np.interp(curr_time, timevalues, resvalues[j])
                obj[j].move_to(np.sin(u) * np.cos(phi(u)) * RIGHT + np.sin(u) * np.sin(phi(u)) * UP + np.cos(u) * OUT)
        
        pdots.add_updater(pdotsupdater)
        self.add(pdots)

        self.play(
            t.animate.set_value(4 * total_time),
            #self.camera.animate.set_theta(390 * DEGREES),
            run_time = total_time,
            rate_func = linear
        )

        pdots.remove_updater(pdotsupdater)
        #self.stop_ambient_camera_rotation()
        self.wait(1)
        self.interactive_embed()


class SphericalTautochroneCairo(ThreeDScene):
    def construct(self):
        #https://sam-dolan.staff.shef.ac.uk/mas212/notebooks/ODE_Example.html
        #https://pastebin.com/2s9At2XE
        sph_rad = 2

        t = ValueTracker(0)
        tlab = DecimalNumber(0)
        tlab.set_value(t.get_value(), num_decimal_places = 2)
        tlab.add_updater(lambda m: m.set_value(t.get_value()))
        tgp = VGroup(MathTex("t="),tlab)
        tgp.arrange(RIGHT)

        total_time = 6
        g = 9.81
        lmda = 2 * g * (total_time / PI) * (total_time / PI)
        lmda = 1 / lmda
        lmda = 1 - lmda
        lmda = np.arccos(lmda)

        def phi(theta):
            ct, cl = np.cos(theta), np.cos(lmda)
            return (1 / np.tan(lmda / 2)) * np.arctanh(np.sqrt((cl - ct) / (1 + cl))) - np.arctan(np.sqrt((cl - ct) / (1 - cl)))
        
        def spherical_diff(theta, ti, start_angle):
            sa = max(start_angle, lmda)
            ct, cs = np.cos(theta), np.cos(sa - 1 / 1024 / 1024)
            return (PI / total_time) * np.sqrt((cs - ct) / (1 - ct))
        
        timevalues = np.linspace(0, total_time, 500)

        sph = Sphere(radius = sph_rad).fade(1 / 2)
        ax = ThreeDAxes()
        base_plane = NumberPlane().shift(sph_rad * IN)
        self.add(sph, base_plane)
        self.move_camera(phi = 60 * DEGREES, theta = 30 * DEGREES)

        sph_tauto = ParametricFunction(
            lambda u: sph_rad * np.array([np.sin(u) * np.cos(phi(u)), np.sin(u) * np.sin(phi(u)), np.cos(u)]), 
            t_range = np.array([lmda, PI - 1 / 1024]), fill_opacity = 0, color = RED
        )
        self.play(Write(sph_tauto))

        n = 3
        pdots = VGroup(*[Dot3D(radius = 1 / 16) for k in range(n)])
        gps = (PI - lmda) / n
        resvalues = []
        for j in range(n):
            u = lmda + j * gps
            vals = odeint(spherical_diff, u, timevalues, args = (u, ))
            tvals = vals[:, 0]
            resvalues += [tvals]
            pdots[j].move_to(sph_rad * (np.sin(u) * np.cos(phi(u)) * RIGHT + np.sin(u) * np.sin(phi(u)) * UP + np.cos(u) * OUT))
        
        self.play(
            Write(pdots),
        )
        self.begin_ambient_camera_rotation(rate = 1 / 256)

        def pdotsupdater(obj):
            curr_time = t.get_value() / 4
            for j in range(n):
                u = np.interp(curr_time, timevalues, resvalues[j])
                obj[j].move_to(sph_rad * (np.sin(u) * np.cos(phi(u)) * RIGHT + np.sin(u) * np.sin(phi(u)) * UP + np.cos(u) * OUT))
        
        pdots.add_updater(pdotsupdater)
        self.add(pdots)

        self.play(
            t.animate.set_value(4 * total_time),
            run_time = total_time,
            rate_func = linear
        )

        pdots.remove_updater(pdotsupdater)
        self.stop_ambient_camera_rotation()
        self.wait(1)
        self.interactive_embed()


class BrachisIntro(MovingCameraScene):
    def construct(self):
        #https://sam-dolan.staff.shef.ac.uk/mas212/notebooks/ODE_Example.html
        #https://pastebin.com/2s9At2XE
        t = ValueTracker(0)
        tlab = DecimalNumber(0)
        tlab.set_value(t.get_value(), num_decimal_places = 2)
        tlab.add_updater(lambda m: m.set_value(t.get_value()))
        tgp = VGroup(MathTex("t="),tlab)
        tgp.arrange(RIGHT)

        rad = 3
        tracing_circle = Circle(radius = rad, color = WHITE)
        path, tracing_dot = VMobject(), Dot(color = RED, stroke_width = 3, stroke_color = WHITE)
        n_sects = 7
        sects = AnnularSector(inner_radius = 0, outer_radius = rad, fill_opacity = 1 / 2, angle = 360 * DEGREES / n_sects)
        tracing_grp = VGroup()
        cols = [RED, ORANGE, GOLD, YELLOW, GREEN, BLUE, PURPLE]
        for k in range(n_sects):
            temp = sects.copy()
            temp.set_style(fill_opacity = 1 / 2, fill_color = cols[k])
            temp.rotate_about_origin(k * 360 * DEGREES / n_sects)
            tracing_grp.add(temp)
        tracing_grp.add(tracing_circle)
        tracing_grp.rotate_about_origin(90 * DEGREES)
        tracing_dot.move_to(tracing_grp[-1].get_start())

        tgrp = VGroup(tracing_grp, tracing_dot)
        path.set_points_as_corners([tracing_dot.get_center(), tracing_dot.get_center()])
        tgrp.save_state()

        j = ValueTracker(0)

        def update_tracing_grp(obj):
            tgrp.restore()
            tgrp.rotate_about_origin(j.get_value())
            tgrp.shift(j.get_value() * rad * RIGHT)
            obj.become(tgrp)
        
        def update_path(path):
            previous_path = path.copy()
            previous_path.add_points_as_corners([tracing_dot.get_center()])
            path.become(previous_path)
        
        path.add_updater(update_path)
        tgrp.add_updater(update_tracing_grp)

        self.camera.frame.shift(5 * RIGHT)

        self.play(Write(tgrp))
        self.add(path, tgrp)

        self.play(
            j.animate.set_value(PI),
            run_time = 5,
            rate_func = linear
        )
        
        num_pts = 1024
        step_size = rad * PI / num_pts
        cycloid_points = [rad * (u - np.sin(u)) * RIGHT + rad * (1 - np.cos(u)) * DOWN + rad * UP for u in np.arange(0, PI + step_size, step_size)]
        cycloid_curve = VMobject()
        cycloid_curve.set_points_as_corners(cycloid_points)
        cycloid_points += [rad * DOWN, rad * UP]
        self.add(cycloid_curve)

        path.remove_updater(update_path)
        tgrp.remove_updater(update_tracing_grp)

        pedestal = VMobject(fill_opacity = 1 / 2, fill_color = GREY)
        pedestal.set_points_as_corners(cycloid_points)

        self.play(FadeOut(tgrp), Write(pedestal))
        self.remove(cycloid_curve)
        self.wait()

        other_line = Line(rad * UP, rad * DOWN + rad * PI * RIGHT)
        line_dot, cycloid_dot = Dot(rad * UP, color = GREEN, stroke_width = 3, stroke_color = WHITE), Dot(rad * UP, color = YELLOW, stroke_width = 3, stroke_color = WHITE)
        self.play(Write(other_line), Write(line_dot), Write(cycloid_dot))
        self.wait()

        cycloid_inits = [1 * DEGREES / 1024, 0.0]
        g, line_incline = 9.81, np.arctan(rad * PI / 2 / rad)
        total_time = np.sqrt(2 * other_line.get_length() / g / np.cos(line_incline))
        timevalues = np.linspace(0, total_time, 1024)
        dots_grp = VGroup(line_dot, cycloid_dot)

        def cycloid_diffeqn(vec, t):
            return [vec[1], (g / 2 / 3 - vec[1] * vec[1] / 2) / np.tan(vec[0] / 2)]
        res = odeint(cycloid_diffeqn, cycloid_inits, timevalues)
        cycloid_values = res[:, 0]
        t.set_value(0)

        def update_dots(obj):
            curr_time = t.get_value() / 4
            u = np.interp(curr_time, timevalues, cycloid_values)
            u = min(u, PI)
            obj[1].move_to(rad * (u - np.sin(u)) * RIGHT + rad * (1 - np.cos(u)) * DOWN + rad * UP)
            obj[0].move_to(g * np.cos(line_incline) * curr_time * curr_time / 2 * DOWN)
            obj[0].rotate_about_origin(line_incline)
            obj[0].shift(rad * UP)
        
        dots_grp.add_updater(update_dots)
        self.add(dots_grp)

        brachis_label = MathTex("\\text{Brachisto}", "\\cdot", "\\text{chrone}").move_to(0.75 * RIGHT * rad * PI + 0.75 * rad * UP)
        brachis_braces = VGroup(
            Brace(brachis_label[0], DOWN),
            Brace(brachis_label[-1], DOWN),
        )
        brachis_meaning = VGroup(Tex("shortest"), Tex("time")).set_color(YELLOW)
        brachis_meaning[0].next_to(brachis_braces[0], DOWN)
        brachis_meaning[1].next_to(brachis_braces[1], DOWN)

        '''self.play(
            AnimationGroup(
                Write(brachis_label),
                Write(VGroup(brachis_braces, brachis_meaning)),
                run_time = 1 / 2
            ),
            t.animate.set_value(4 * total_time),
            run_time = 4 * total_time,
            rate_func = linear
        )'''

        self.play(
            AnimationGroup(
                ApplyMethod(t.set_value, 4 * total_time, run_time = 4 * total_time, rate_func = linear),
                AnimationGroup(Write(brachis_label), Write(VGroup(brachis_braces, brachis_meaning)), run_time = 2, lag_ratio = 1),
                lag_ratio = 1 / 4
            )
        )

        dots_grp.remove_updater(update_dots)
        self.play(FadeOut(dots_grp), FadeOut(other_line))
        self.wait()

        seven_dots = VGroup()
        total_time = PI * np.sqrt(rad / g)
        timevalues = np.linspace(0, total_time, 1024)
        seven_dots_values = []
        for k in range(n_sects):
            temp = Dot(rad * UP, color = cols[k], stroke_width = 3, stroke_color = WHITE)
            u = k * PI / n_sects
            temp.move_to(rad * (u - np.sin(u)) * RIGHT + rad * (1 - np.cos(u)) * DOWN + rad * UP)
            seven_dots.add(temp)
            temp_init = [max(u, 1 * DEGREES / 1024), 0]
            res = odeint(cycloid_diffeqn, temp_init, timevalues)
            seven_dots_values += [res[:, 0]]
        
        tauto_label = Tex("Tauto").next_to(brachis_label[1], LEFT)
        tauto_label.shift(UP / 16 + DOWN / 32 + RIGHT / 16)
        tauto_brace = Brace(tauto_label, DOWN)
        tauto_meaning = Tex("same").set_color(YELLOW).next_to(tauto_brace, DOWN)
        
        self.play(
            Write(seven_dots),
            Transform(brachis_label[0], tauto_label),
            Transform(brachis_braces[0], tauto_brace),
            Transform(brachis_meaning[0], tauto_meaning)
        )
        self.wait()

        t.set_value(0)
        def update_seven_dots(obj):
            curr_time = t.get_value() / 4
            for j in range(n_sects):
                u = min(np.interp(curr_time, timevalues, seven_dots_values[j]), PI)
                obj[j].move_to(rad * (u - np.sin(u)) * RIGHT + rad * (1 - np.cos(u)) * DOWN + rad * UP)
        seven_dots.add_updater(update_seven_dots)
        self.add(seven_dots)

        self.play(
            t.animate.set_value(4 * total_time),
            run_time = 4 * total_time,
            rate_func = linear
        )
        seven_dots.remove_updater(update_seven_dots)
        self.wait()

        self.wait(5)


class SimplePendulum(MovingCameraScene):
    def construct(self):
        #https://sam-dolan.staff.shef.ac.uk/mas212/notebooks/ODE_Example.html
        #https://pastebin.com/2s9At2XE
        t = ValueTracker(0)
        tlab = DecimalNumber(0)
        tlab.set_value(t.get_value(), num_decimal_places = 2)
        tlab.add_updater(lambda m: m.set_value(t.get_value()))
        tgp = VGroup(MathTex("t="),tlab)
        tgp.arrange(RIGHT)

        self.camera.frame.shift(3 * DOWN)

        l, rad, g = 6, 1 / 2, 9.81
        omeg = np.sqrt(g / l)
        #blob = Dot(radius = 1, color = color_gradient([GREY, WHITE, GRAY], 3), stroke_width = 3, stroke_color = WHITE)
        blob = Dot(radius = rad, fill_color = color_gradient([GRAY, WHITE], 4))
        pendulum_string = Line(ORIGIN, l * DOWN)
        blob.move_to(pendulum_string.get_end())
        pend = VGroup(pendulum_string, blob)

        wall = Line(LEFT / 2, RIGHT / 2)
        wline = Line(ORIGIN, RIGHT / 4, color = GREY).rotate_about_origin(45 * DEGREES)
        nwlines = 8
        wall_gp = VGroup()
        pendulum_pivot = Dot()
        for k in range(nwlines):
            temp = wline.copy().shift(k * RIGHT / nwlines + LEFT / 2)
            wall_gp.add(temp)
        wall_gp.add(pendulum_pivot, wall)

        self.play(Write(pend), Write(wall_gp))

        velocity_init = 1 / 2

        def update_pendulum(obj):
            tp, tb = obj
            ctime = t.get_value() / 4
            tp.become(Line(ORIGIN, l * DOWN))
            ctheta = velocity_init * np.sin(omeg * ctime) / omeg
            tp.rotate_about_origin(ctheta)
            tb.move_to(tp.get_end())
        pend.add_updater(update_pendulum)
        self.add(pend)

        time_period = 2 * PI / omeg
        cycles = 3

        self.play(
            t.animate.set_value(4 * cycles * time_period),
            run_time = cycles * time_period,
            rate_func = linear
        )
        self.wait()

        pend.remove_updater(update_pendulum)
        t.set_value(0)

        self.play(self.camera.frame.animate.shift(4 * RIGHT))
        self.wait()

        period_label = MathTex("\\text{period of motion, }", "T_0", "=", "\\frac{\\pi}{2}", "\\sqrt{\\frac{l}{g}}")
        period_label[1].set_color(BLUE)
        period_label[-1][2].set_color(ORANGE)
        period_label[-1][-1].set_color(YELLOW)
        period_label.move_to(6 * RIGHT + 1.5 * DOWN)

        self.play(Write(period_label))
        self.wait()

        pend.add_updater(update_pendulum)
        self.play(
            t.animate.set_value(-4 * PI / 2 / omeg),
            run_time = 2,
            rate_func = linear
        )
        self.wait()

        self.play(
            t.animate.set_value(0),
            run_time = 2,
            rate_func = linear
        )
        self.wait()

        pe_label = MathTex("\\text{Potential Energy, }", "V", "=", "m", "g", "l", "(1-\\cos\\theta)")
        pe_label[-2].set_color(ORANGE)
        pe_label[-3].set_color(YELLOW)
        pe_label[-1][-2].set_color(GREEN)
        pe_label.next_to(period_label, DOWN, buff = 1)
        pe_label.shift(1.5 * RIGHT)
        misc_objects = VGroup(
            Line(blob.get_center() + 3 * LEFT, blob.get_center() + 3 * RIGHT).set_color(GREY),
            Line(ORIGIN, l * DOWN).set_color(GREY),
        )

        self.play(
            t.animate.set_value(1.0 * 4 * PI / 2 / omeg),
            Write(misc_objects)
        )
        self.wait()
        pend.remove_updater(update_pendulum)

        ang = Anglet(misc_objects[1].get_end(), ORIGIN, pendulum_string.get_end(), radius = 1)
        theta_label = MathTex("\\theta").set_color(GREEN)
        theta_label.next_to(ang, RIGHT, buff = 1 / 8)
        hline = Line(blob.get_center(), blob.get_center() + 3 * LEFT).set_color(GREY)
        curr_y = blob.get_center()[1]
        lcosbrace = Brace(Line(ORIGIN, curr_y * UP), LEFT)
        lcos_label = MathTex("l\\cos\\theta").rotate_in_place(90 * DEGREES)
        lcos_label[0][-1].set_color(GREEN)
        lcos_label.next_to(lcosbrace, LEFT)

        self.play(
            self.camera.frame.animate.shift(1.5 * RIGHT),
            period_label.animate.shift(1.5 * RIGHT),
            Write(ang),
            Write(theta_label),
            Write(pe_label),
            Write(hline),
            Write(lcosbrace),
            Write(lcos_label)
        )
        self.wait()

        small_cos = MathTex("\\frac{\\theta^2}{2}")
        small_cos[0][0].set_color(GREEN)
        small_cos.move_to(pe_label[5].get_right() + small_cos.get_center() - small_cos[0][0].get_left() + RIGHT / 16)
        #small_cos.next_to(pe_label[5], RIGHT)

        self.play(
            FadeOut(VGroup(lcos_label, lcosbrace, hline)),
            Transform(pe_label[-1], small_cos)
        )
        self.wait()

        alternate_pe = MathTex("\\frac{1}{2}", "m", "\\frac{g}{l}", "(l\\theta)", "^2")
        alternate_pe[2][0].set_color(YELLOW)
        alternate_pe[2][-1].set_color(ORANGE)
        alternate_pe[-2][1].set_color(ORANGE)
        alternate_pe[-2][2].set_color(GREEN)
        alternate_pe.next_to(pe_label[2], RIGHT)

        self.play(ReplacementTransform(pe_label[3:], alternate_pe))
        self.wait()

        length_arc = Arc(radius = l, angle = -velocity_init * np.sin(4 * PI / 2 / omeg) / omeg).set_color(GREY)
        length_brace = ArcBrace(length_arc)
        length_gp = VGroup(length_arc, length_brace)
        length_gp.rotate_about_origin(-90 * DEGREES)

        slabel = MathTex("s").set_color(GREEN)
        slabel.next_to(length_brace, DOWN)

        self.play(
            FadeOut(VGroup(ang, theta_label, misc_objects[0])),
            FadeIn(length_gp),
            Write(slabel)
        )
        self.wait()

        pe_slabel = slabel.copy()
        pe_slabel.move_to(alternate_pe[-2].get_center())

        self.play(
            Transform(alternate_pe[-2], pe_slabel)
        )
        self.wait()

        gbyl = MathTex("\\frac{g}{l}=\\left(\\frac{\\pi}{2T_0}\\right)^2")
        gbyl[0][0].set_color(YELLOW)
        gbyl[0][2].set_color(ORANGE)
        gbyl[0][8:10].set_color(BLUE)
        gbyl.move_to(period_label.get_center())

        self.play(
            FadeOut(period_label, shift = UP),
            FadeIn(gbyl, shift = UP)
        )
        self.wait()

        new_pe = MathTex("V", "=", "\\frac{1}{2}", "m", "\\left(\\frac{\\pi}{2T_0}\\right)^2", "s^2")
        new_pe[-1][-2].set_color(GREEN)
        new_pe[-2][4:6].set_color(BLUE)

        self.play(
            FadeOut(VGroup(gbyl, length_gp, slabel, pend, wall_gp, misc_objects[1], pe_label)),
            self.camera.frame.animate.move_to(ORIGIN),
            Write(new_pe)
        )
        self.wait()

        arb_cur = VGroup(
            Line(DOWN / 4, 3 * UP).set_color(GREY),
            Line(LEFT / 4, 4 * RIGHT).set_color(GREY),
            ArcBetweenPoints(1 * UR, 1 * UR + 3 * UR).set_color(BLUE)
        )
        rpt = Dot(arb_cur[2].point_from_proportion(0.5), color = RED)
        ax_labs = VGroup(MathTex("x"), MathTex("y"))
        ax_labs[0].next_to(arb_cur[1].get_end(), RIGHT)
        ax_labs[1].next_to(arb_cur[0].get_end(), UP)
        arb_cur.add(rpt)
        arb_cur.add(ax_labs)
        arb_cur.shift(2.5 * LEFT)

        yline = Line(rpt.get_center(), rpt.get_center()[0] * RIGHT, color = GREY)
        ybrace = Brace(yline, RIGHT)
        ylabel = MathTex("y").next_to(ybrace, RIGHT).set_color(ORANGE)
        arb_cur.add(yline, ybrace, ylabel)

        self.play(
            new_pe.animate.shift(1.5 * DOWN),
            Write(arb_cur[: -3]),
        )
        self.wait()

        newpec = new_pe.copy()
        newpec.shift(LEFT)
        pe_rhs = MathTex("=", "m", "g", "y")
        pe_rhs.next_to(newpec, RIGHT).shift(DOWN / 8)
        pe_rhs[-1].set_color(ORANGE)
        pe_rhs[-2].set_color(YELLOW)

        self.play(
            new_pe.animate.shift(LEFT),
            FadeIn(pe_rhs, shift = LEFT),
            Write(arb_cur[-3:])
        )
        self.wait()

        sequals = MathTex("s", "=", "\\pm", "\\frac{2T_0}{\\pi}", "\\sqrt{2gy}")
        sequals.to_edge(DOWN)
        sequals[0].set_color(GREEN)
        sequals[3][1:3].set_color(BLUE)
        sequals[-1][-1].set_color(ORANGE)
        sequals[-1][-2].set_color(YELLOW)

        self.play(Write(sequals))
        self.wait()

        self.play(
            FadeOut(VGroup(arb_cur, new_pe, pe_rhs), shfit = UP),
            sequals.animate.move_to(2 * UP)
        )
        self.wait()

        s_diff = MathTex("\\frac{ds}{dy}", "=", "-", "\\frac{T_0}{\\pi}", "\\sqrt{\\frac{2g}{y}}")
        s_diff[0][1].set_color(GREEN)
        s_diff[0][-1].set_color(ORANGE)
        s_diff[3][0:2].set_color(BLUE)
        s_diff[-1][-1].set_color(ORANGE)
        s_diff[-1][-3].set_color(YELLOW)

        self.play(Write(s_diff))
        self.wait()

        surr = SurroundingRectangle(s_diff)
        self.play(Create(surr))

        difftauto = Tex("Differential Equation of a ", "Tautochrone", " Curve")
        difftauto[-2].set_color(BLUE)
        diffuline = Underline(difftauto)
        difftgp = VGroup(difftauto, diffuline)
        difftgp.to_edge(UP)

        self.play(
            FadeOut(sequals),
            Write(difftgp)
        )

        self.wait(5)


class SphericalCoordinates(ThreeDScene):
    def construct(self):
        #https://sam-dolan.staff.shef.ac.uk/mas212/notebooks/ODE_Example.html
        #https://pastebin.com/2s9At2XE
        sph_rad = 2.5
        sph = Sphere(radius = sph_rad)
        ax = ThreeDAxes()
        base_plane = NumberPlane().shift(sph_rad * IN)
        xyplane = VMobject(fill_color = BLUE, stroke_width = 1 / 4, fill_opacity = 7 / 8)
        xyplane.set_points_as_corners([5 * UL, 5 * UR, 5 * DR, 5 * DL, 5 * UL])
        self.add(ax)

        self.play(Write(sph))
        self.move_camera(phi = 60 * DEGREES, theta = 30 * DEGREES)
        self.wait()

        #self.play(sph.animate.set_style(fill_opacity = 1 / 8))

        pdot = Dot3D(sph_rad * OUT)
        pline = Line(ORIGIN, pdot.get_center())
        self.play(Write(pdot), Write(pline))
        self.wait()

        polar_grid = VGroup()
        circle_mesh = VGroup()
        nums = 4
        for k in range(1, nums + 1):
            rd = k * sph_rad / nums
            a = k * 90 * DEGREES / nums
            circle_mesh.add(
                Arc(radius = rd, stroke_width = 1 / 2),
                Line(ORIGIN, sph_rad * RIGHT, stroke_width = 1 / 2).rotate_about_origin(a)
            )
        #circle_mesh.set_color(GREY)
        polar_grid.add(circle_mesh.copy())
        polar_grid.add(circle_mesh.copy().rotate(90 * DEGREES, about_point = ORIGIN, axis = RIGHT))
        polar_grid.add(circle_mesh.copy().rotate(-90 * DEGREES, about_point = ORIGIN, axis = UP))
        #carc = Sector(outer_radius = sph_rad, fill_color = GREEN, fill_opacity = 1 / 2, stroke_width = 2)
        carc = VMobject(fill_color = BLUE, fill_opacity = 1 / 4, stroke_color = WHITE)
        cpts = [ORIGIN]
        for k in range(100 + 1):
            cpts += [sph_rad * np.cos(k * PI / 2 / 100) * RIGHT + sph_rad * np.sin(k * PI / 2 / 100) * UP]
        cpts += [ORIGIN]
        carc.set_points_as_corners(cpts)
        polar_grid.add(carc.copy())
        polar_grid.add(carc.copy().rotate(90 * DEGREES, about_point = ORIGIN, axis = RIGHT))
        polar_grid.add(carc.copy().rotate(-90 * DEGREES, about_point = ORIGIN, axis = UP))

        self.play(
            sph.animate.fade(7 / 8),
            FadeIn(polar_grid)
        )
        #self.move_camera(theta = 90 * DEGREES)
        self.wait()

        ta, pa = ValueTracker(0), ValueTracker(0)
        #pdot.add_updater(lambda m: m.restore().rotate().rotate())
        pdot.add_updater(lambda m: m.move_to(sph_rad * (np.cos(ta.get_value()) * OUT + np.sin(ta.get_value()) * np.cos(pa.get_value()) * RIGHT + np.sin(ta.get_value()) * np.sin(pa.get_value()) * UP)))
        pline.add_updater(lambda m: m.become(Line(ORIGIN, pdot.get_center())))
        self.add(pdot, pline)

        self.play(ta.animate.set_value(45 * DEGREES), run_time = 2)
        self.play(pa.animate.set_value(45 * DEGREES), run_time = 2)
        self.wait()

        xydot, xdot = Dot(sph_rad * np.sin(ta.get_value()) * (np.cos(pa.get_value()) * RIGHT + np.sin(pa.get_value()) * UP)), Dot(sph_rad * np.sin(ta.get_value() * np.cos(pa.get_value()) * RIGHT), fill_opacity = 0)
        zline = Line(pdot.get_center(), xydot.get_center(), color = GREY)
        yline = Line(xydot.get_center(), xdot.get_center(), color = GRAY)
        xyline = Line(ORIGIN, xydot.get_center(), color = GREY)
        
        self.play(Write(zline), Write(yline), Write(xyline), Write(xydot), Write(xdot))
        self.wait()

        xydot.add_updater(lambda m: m.move_to(sph_rad * np.sin(ta.get_value()) * (np.cos(pa.get_value()) * RIGHT + np.sin(pa.get_value()) * UP)))
        xdot.add_updater(lambda m: m.move_to(sph_rad * np.sin(ta.get_value() * np.cos(pa.get_value()) * RIGHT)))
        zline.add_updater(lambda m: m.become(Line(pdot.get_center(), xydot.get_center(), color = GREY)))
        yline.add_updater(lambda m: m.become(Line(xydot.get_center(), xdot.get_center(), color = GREY)))
        xyline.add_updater(lambda m: m.become(Line(ORIGIN, xydot.get_center(), color = GREY)))
        self.add(xydot, xdot, zline, yline, xyline)

        '''thet_ang, phi_ang = Anglet(pdot.get_center(), ORIGIN, OUT, color = GREEN), Anglet(xydot.get_center(), ORIGIN, xdot.get_center(), color = RED)
        self.play(Write(thet_ang), Write(phi_ang))
        thet_ang.add_updater(lambda m: m.become(Anglet(pdot.get_center(), ORIGIN, OUT, color = GREEN)))
        phi_ang.add_updater(lambda m: m.become(Anglet(xydot.get_center(), ORIGIN, xdot.get_center(), color = RED)))
        self.add(thet_ang, phi_ang)'''

        thet_circle = Circle(color = GREEN, radius = sph_rad).rotate_about_origin(PI / 2, axis = RIGHT).rotate_about_origin(pa.get_value())
        phi_circle = Circle(color = RED, radius = sph_rad * np.sin(ta.get_value())).shift(sph_rad * np.cos(ta.get_value()) * OUT)

        angle_defs = VGroup(
            MathTex("\\text{polar angle, }\\theta"),
            MathTex("\\text{azimuthal angle, }\\phi"),
        )
        angle_defs.arrange(DOWN)
        angle_defs[0][0][-1].set_color(GREEN)
        angle_defs[1][0][-1].set_color(RED)
        angle_defs.to_corner(UR)

        self.add_fixed_in_frame_mobjects(angle_defs[0])
        self.play(Write(thet_circle))
        self.wait()

        thet_circle.add_updater(lambda m: m.become(Circle(color = GREEN, radius = sph_rad).rotate_about_origin(PI / 2, axis = RIGHT).rotate_about_origin(pa.get_value())))
        self.add(thet_circle)

        self.move_camera(theta = -45 * DEGREES, phi = 90 * DEGREES)
        self.play(
            ta.animate.set_value(179 * DEGREES),
            run_time = 5, rate_func = there_and_back,
        )
        self.move_camera(phi=75 * DEGREES, theta=30 * DEGREES)
        self.wait()

        self.add_fixed_in_frame_mobjects(angle_defs[1])
        self.play(Write(phi_circle))
        phi_circle.add_updater(lambda m: m.become(Circle(color = RED, radius = sph_rad * np.sin(ta.get_value())).shift(sph_rad * np.cos(ta.get_value()) * OUT)))
        self.add(phi_circle)
        self.wait()

        self.move_camera(phi = 0 * DEGREES)
        self.play(pa.animate.increment_value(180 * DEGREES), run_time = 5, rate_func = there_and_back)
        self.move_camera(phi = 45 * DEGREES)
        self.wait()

        z_arclen = VGroup(
            MathTex("z = 1 + \\cos \\theta"),
            MathTex("ds^2 = d\\theta^2 + \\sin^2\\theta d\\phi^2")
        )
        z_arclen.arrange(DOWN)
        z_arclen[0][0][0].set_color(ORANGE)
        z_arclen[0][0][-1].set_color(GREEN)
        z_arclen[1][0][1].set_color(GOLD)
        z_arclen[1][0][5].set_color(GREEN)
        z_arclen[1][0][12].set_color(GREEN)
        z_arclen[1][0][14].set_color(RED)
        z_arclen.add_background_rectangle()

        z_arclen.to_corner(DR)

        linegp = VGroup(
            Line(sph_rad * DOWN, sph_rad * np.cos(ta.get_value()) * UP, color = ORANGE),
        )
        linegp.add(Brace(linegp[-1], LEFT))
        linegp.add(MathTex("1+\\cos\\theta").rotate_about_origin(90 * DEGREES).next_to(linegp[-1], LEFT))
        linegp[-1][0][-1].set_color(GREEN)
        linegp.add(Line(3 * LEFT + sph_rad * DOWN, 3 * RIGHT + sph_rad * DOWN))
        linegp.rotate_about_origin(90 * DEGREES, axis = RIGHT)
        linegp.rotate_about_origin(pa.get_value())
        linegp.add(Line(pdot.get_center(), linegp[0].get_end(), color = GREY))

        self.move_camera(theta = -PI / 2 + pa.get_value(), phi = 90 * DEGREES, run_time = 2)
        self.play(Write(linegp))
        self.add_fixed_in_frame_mobjects(z_arclen)
        self.wait()

        self.move_camera(theta = 15 * DEGREES, phi = 60 * DEGREES, run_time = 2)

        self.begin_ambient_camera_rotation(rate = 0.1)

        xydot.clear_updaters()
        xdot.clear_updaters()
        pdot.clear_updaters()
        pline.clear_updaters()
        zline.clear_updaters()
        yline.clear_updaters()
        xyline.clear_updaters()
        #thet_ang.clear_updaters()
        #phi_ang.clear_updaters()
        thet_circle.clear_updaters()
        phi_circle.clear_updaters()
        self.wait(10)


class FinalDiffEqn(MovingCameraScene):
    def construct(self):
        diffeqn = MathTex("ds", "=", "\\sqrt{2g}", "{T_0", "\\over", "\\pi}", "{dz", "\\over", "\\sqrt{z}}")
        diffeqn[0][-1].set_color(GREEN)
        diffeqn[2][-1].set_color(YELLOW)
        diffeqn[3].set_color(BLUE)
        diffeqn[6][-1].set_color(ORANGE)
        diffeqn[-1][-1].set_color(ORANGE)
        self.play(Write(diffeqn))
        self.wait()
        self.play(diffeqn.animate.to_edge(UP))
        self.wait()

        vals = VGroup(
            MathTex("ds^2", "=", "d\\theta^2", "+", "\\sin\\theta^2d\\phi"),
            MathTex("2g", "\\left(\\frac{T_0}{\\pi}\\right)^2", "=", "\\frac{1}{1-\\cos\\lambda}"),
            MathTex("z", "=", "1+\\cos\\theta"),
        )
        vals[0][0][1].set_color(GREEN)
        vals[0][2][1].set_color(ORANGE)
        vals[0][4][3].set_color(ORANGE)
        vals[0][4][-1].set_color(TEAL)
        vals[1][0][-1].set_color(YELLOW)
        vals[1][1][1:3].set_color(BLUE)
        vals[1][3][-1].set_color(BLUE)
        vals[2][0][0].set_color(ORANGE)
        vals[2][2][-1].set_color(ORANGE)

        vals[0].shift(3 * LEFT )
        vals[2].shift(3 * RIGHT)
        vals[1].shift(1.5 * DOWN)

        self.play(
            Write(vals[0]),
            Write(vals[2]),
        )
        self.wait()
        
        self.play(
            vals[0].animate.shift(UP / 2),
            vals[2].animate.shift(UP / 2),
            Write(vals[1])
        )
        self.wait()

        phidiff = MathTex("{d\\phi", "\\over", "d\\theta}", "=", "{1", "\\over", "\\sin\\theta}", "\\sqrt{", "{{\\cos\\lambda", "-", "\\cos\\theta}", "\\over", "{1-\\cos\\lambda}}}")
        phidiff[0][-1].set_color(TEAL)
        phidiff[2][-1].set_color(ORANGE)
        phidiff[6][-1].set_color(ORANGE)
        phidiff[8][-1].set_color(BLUE)
        phidiff[10][-1].set_color(ORANGE)
        phidiff[12][-1].set_color(BLUE)

        self.play(
            FadeOut(vals[0], shift = LEFT),
            FadeOut(vals[2], shift = RIGHT),
            FadeOut(vals[1], shift = DOWN),
            FadeOut(diffeqn, shift = UP),
            Write(phidiff)
        )
        self.wait()

        ques2 = Tex("Looks very hard..")
        ques1 = Tex(".. But can be solved \\\\ in closed form")
        bub = SVGMobject("/Users/muthuveerappanramalingam/Downloads/ManimInstall/manim_ce/GeometricProbability/assets/Bubbles_thought")
        bub.set_color(WHITE).set_style(fill_color = BLACK, fill_opacity = 3 / 4, stroke_color = WHITE, stroke_width = 3).scale_in_place(2).stretch_in_place(2, 0)
        bub.flip()
        def add_and_resize(bubb, contt, content_scale_factor = 0.75, bub_cent_adj_fac = 1 / 8):
            bub, cont = bubb.copy(), contt
            twid = cont.width
            twid += max(2, MED_LARGE_BUFF)
            tht = cont.height
            tht += 2.5 * LARGE_BUFF
            bub.width, bub.height = twid, tht
            bubcent = bub.get_center() + bub_cent_adj_fac * bub.height * UP
            swid = content_scale_factor * bub.width
            if cont.width > swid:
                cont.width = swid
            cont.shift(bubcent - cont.get_center())
            return VGroup(bub, cont)
        ques1grp = add_and_resize(bub.copy(), ques1)
        ques2grp = add_and_resize(bub.copy().flip(), ques2)
        ques1grp.to_corner(DR)
        ques2grp.to_corner(DL)

        self.play(
            phidiff.animate.shift(3 * UP),
            Write(ques2grp)
        )
        self.play(Write(ques1grp))
        self.wait()

        finalphi = MathTex(
            "\\phi", "=",
            "\\cot\\left(\\frac{\\lambda}{2}\\right)",
            "\\tanh^{-1}", "\\left(", "\\sqrt{", "{\\cos\\lambda", "-", "\\cos\\theta}", "\\over", "{1-\\cos\\lambda}}", "\\right)",
            "-",
            "\\tan^{-1}", "\\left(", "\\sqrt{", "{\\cos\\lambda", "-", "\\cos\\theta}", "\\over", "{1+\\cos\\lambda}}", "\\right)"
        )
        finalphi.scale(15 / 16)
        finalphi[0].set_color(TEAL)
        finalphi[2][-4].set_color(BLUE)
        finalphi[6][-1].set_color(BLUE)
        finalphi[8][-1].set_color(ORANGE)
        finalphi[10][-1].set_color(BLUE)
        finalphi[16][-1].set_color(BLUE)
        finalphi[18][-1].set_color(ORANGE)
        finalphi[20][-1].set_color(BLUE)

        title = Tex("Tautochrone", " on a ", "Sphere")
        title.shift(2 * UP)
        titleu = Underline(title)
        titlegp = VGroup(title, titleu)
        title[0].set_color(GOLD)
        title[-1].set_color(YELLOW)

        self.play(
            FadeOut(ques1grp, shift = LEFT),
            FadeOut(ques2grp, shift = RIGHT),
            Write(finalphi),
            FadeOut(phidiff, shift = UP)
        )
        self.play(Write(titlegp))

        self.wait(5)


if __name__ == "__main__":
    module_name = os.path.abspath(__file__)
    #output_location = "C:\ManimCE\media"
    #clear_cmd = "cls"
    #command_A = "manim " + module_name + " " + "RealCase" + " " + "-pql -n 42" + " --media_dir " + output_location
    #command_A = "manim " + module_name + " --media_dir " + output_location + " -pqh"
    #command_A = "manim "+ "-p" + " " + module_name + " " + "SphericalTautochrone" + " -n 0"+ " --renderer=opengl"
    #command_A = "manim "+ "-pql" + " " + module_name + " " + "SphericalTautochroneCairo" + " -n 0"
    #command_A = "manim "+ "-pqh" + " " + module_name
    command_A = "manim "+ "--write_to_movie" + " " + module_name + " " + "SphericalTautochrone" + " -n 0"+ " --renderer=opengl"
    #os.system(clear_cmd)
    os.system(command_A)
