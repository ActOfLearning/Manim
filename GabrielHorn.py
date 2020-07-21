from manimlib.imports import *
import numpy as np


class RoseCurve(ParametricFunction):
    CONFIG = {
        "radius"        : 1,
        "start_theta"   : 0,
        "kval"          : 3 / 2,
        "end_theta"     : 2 * PI,
        "step_size"     : 0.001,
        "density"       : 50 * DEFAULT_POINT_DENSITY_1D,
        "color"         : BLUE,
    }
    def __init__(self, **kwargs):
        digest_config(self, kwargs)
        ParametricFunction.__init__(self, self.pos_func, **kwargs)

    def pos_func(self, t):
        T = self.start_theta + t * (self.end_theta - self.start_theta)
        return self.radius * np.array([
            np.cos(self.kval * T) * np.cos(T),
            np.cos(self.kval * T) * np.sin(T),
            0
        ])

class GabHorn(ParametricSurface):
    CONFIG = {
        "resolution": (12, 24),
        "u_min": 1,
        "u_max": 10,
        "v_min": 0,
        "v_max": 2 * PI,
        "checkerboard_colors": [BLUE, GREEN, YELLOW, ORANGE],
    }

    def __init__(self, **kwargs):
        ParametricSurface.__init__(
            self, self.func, **kwargs
        )

    def func(self, u, v):
        return np.array([
            u,
            np.cos(v) / u,
            np.sin(v) / u
        ])

class FlatPlane(ParametricSurface):
    CONFIG = {
        "resolution": (12, 24),
        "u_min": -1,
        "u_max": 1,
        "v_min": -1,
        "v_max": 1,
        "checkerboard_colors": [BLUE, GREEN, YELLOW, ORANGE],
    }

    def __init__(self, **kwargs):
        ParametricSurface.__init__(
            self, self.func, **kwargs
        )

    def func(self, u, v):
        return np.array([
            u,
            v,
            0
        ])

class HemiSphere(ParametricSurface):
    CONFIG = {
        "resolution": (12, 24),
        "radius": 1,
        "u_min": 0,
        "u_max": PI / 2,
        "v_min": 0,
        "v_max": 2 * PI,
        "checkerboard_colors": [BLUE],
    }

    def __init__(self, **kwargs):
        ParametricSurface.__init__(
            self, self.func, **kwargs
        )
        self.scale(self.radius)

    def func(self, u, v):
        return np.array([
            np.cos(u) * np.cos(v),
            np.cos(u) * np.sin(v),
            np.sin(u)
        ])

class Cone(ParametricSurface):
    CONFIG = {
        "resolution": (12, 24),
        "u_min": 0,
        "u_max": 1,
        "v_min": 0,
        "v_max": TAU,
        "checkerboard_colors": [BLUE_E],
    }

    def __init__(self, **kwargs):
        ParametricSurface.__init__(
            self, self.func, **kwargs
        )
        #self.scale(self.radius)

    def func(self, u, v):
        return np.array([
            np.cos(v) * u,
            np.sin(v) * u,
            u
        ])

class Cylinder(ParametricSurface):
    CONFIG = {
        "resolution": (24, 48),
        "u_min": 0,
        "u_max": 1,
        "v_min": 0,
        "v_max": TAU,
        "radius": 1,
        "checkerboard_colors": [BLUE_E],
    }

    def __init__(self, **kwargs):
        ParametricSurface.__init__(
            self, self.func, **kwargs
        )
        #self.scale(self.radius)

    def func(self, u, v):
        return np.array([
            self.radius * np.cos(v),
            self.radius * np.sin(v),
            u
        ])

class ParametricDisc(Sphere):
    CONFIG = {
        "u_min": 0,
        "u_max": 1,
        "stroke_width": 0,
        "checkerboard_colors": [BLUE_D, BLUE_E],
    }

    def func(self, u, v):
        return np.array([
            u * np.cos(v),
            u * np.sin(v),
            0,
        ])

class Revolution(SpecialThreeDScene):
    def construct(self):
        axes = self.get_axes()
        h = ValueTracker(0.75)
        ang = ValueTracker(1.5 * PI)
        cne = Cone(v_max = ang.get_value(), u_max = h.get_value(), checkerboard_colors = [RED_E])
        cyl = Cylinder(v_max = ang.get_value(), u_max = h.get_value())
        bcone = ParametricDisc(v_max = ang.get_value())
        tcone = ParametricDisc(u_min = h.get_value(), v_max = ang.get_value(), checkerboard_colors = [GREEN_D, GREEN_E])
        tcone.shift(h.get_value() * OUT)
        rwframe = VMobject(fill_color = BLUE, fill_opacity = 0.875)
        rwframe.set_points_as_corners([ORIGIN, RIGHT, RIGHT + h.get_value() * OUT, h.get_value() * RIGHT + h.get_value() * OUT, ORIGIN])
        lwframe = rwframe.copy().rotate_about_origin(ang.get_value(), axis = OUT)
        coneframe = VGroup(
            rwframe, lwframe,
            ParametricFunction(lambda t: np.array([np.cos(t), np.sin(t), h.get_value()]), t_min = 0, t_max = ang.get_value()),
            ParametricFunction(lambda t: np.array([h.get_value() * np.cos(t), h.get_value() * np.sin(t), h.get_value()]), t_min = 0, t_max = ang.get_value()),
        ).set_color(BLUE)
        punccone = VGroup(cne, cyl, bcone, tcone, coneframe)
        punccone.scale(2)
        def cneupd(obj):
            tcne, tcyl, tbcone, ttcone, tconeframe = obj
            tcne.become(Cone(v_max = ang.get_value(), u_max = h.get_value(), checkerboard_colors = [RED_E]))
            tcyl.become(Cylinder(v_max = ang.get_value(), u_max = h.get_value()))
            tbcone.become(ParametricDisc(v_max = ang.get_value()))
            ttcone.become(ParametricDisc(u_min = h.get_value(), v_max = ang.get_value(), checkerboard_colors = [GREEN_D, GREEN_E]))
            ttcone.shift(h.get_value() * OUT)
            trwframe = VMobject(fill_color = BLUE, fill_opacity = 0.875)
            trwframe.set_points_as_corners([ORIGIN, RIGHT, RIGHT + h.get_value() * OUT, h.get_value() * RIGHT + h.get_value() * OUT, ORIGIN])
            tlwframe = trwframe.copy().rotate_about_origin(ang.get_value(), axis = OUT)
            tconeframe.become(
                VGroup(
                    trwframe, tlwframe,
                    ParametricFunction(lambda t: np.array([np.cos(t), np.sin(t), h.get_value()]), t_min = 0, t_max = ang.get_value()),
                    ParametricFunction(lambda t: np.array([h.get_value() * np.cos(t), h.get_value() * np.sin(t), h.get_value()]), t_min = 0, t_max = ang.get_value()),
                ).set_color(BLUE)
            )
            tpunccone = VGroup(tcne, tcyl, tbcone, ttcone, tconeframe)
            tpunccone.scale(2)
            tpunccone.shift(3 * LEFT)
            obj.become(tpunccone)
        punccone.add_updater(cneupd)
        #punccone.set_fill(GREEN_E, opacity = 0.8)
        #punccone.shift(3 * LEFT)
        hem = HemiSphere(u_max = np.arcsin(h.get_value()), v_max = ang.get_value())
        bhem = ParametricDisc(v_max = ang.get_value())
        them = ParametricDisc(u_max = (1 - h.get_value() ** 2) ** 0.5, v_max = ang.get_value(), checkerboard_colors = [GREEN_D, GREEN_E])
        rsframe = VMobject(fill_color = BLUE, fill_opacity = 0.875)
        arcpts = [ORIGIN] + list(Arc(angle = np.arcsin(h.get_value())).points) + [h.get_value() * UP, ORIGIN]
        rsframe.set_points_as_corners(arcpts)
        rsframe.rotate_about_origin(PI / 2, axis = RIGHT)
        lsframe = rsframe.copy().rotate_about_origin(ang.get_value(), axis = OUT)
        hemframe = VGroup(
            lsframe, rsframe,
            ParametricFunction(lambda t: (1 - h.get_value() ** 2) ** 0.5 * np.array([np.cos(t), np.sin(t), 0]) + h.get_value() * OUT, t_min = 0, t_max = ang.get_value()),
        ).set_color(BLUE)
        them.shift(h.get_value() * OUT)
        vhem = VGroup(hem, bhem, them, hemframe)
        vhem.scale(2)
        def hemupd(obj):
            temphem, tempbhem, tempthem, temphemframe = obj
            temphem.become(HemiSphere(u_max = np.arcsin(h.get_value()), v_max = ang.get_value()))
            tempbhem.become(ParametricDisc(v_max = ang.get_value()))
            tempthem.become(ParametricDisc(u_max = (1 - h.get_value() ** 2) ** 0.5, v_max = ang.get_value(), checkerboard_colors = [GREEN_D, GREEN_E]))
            trsframe = VMobject(fill_color = BLUE, fill_opacity = 0.875)
            tarcpts = [ORIGIN] + list(Arc(angle = np.arcsin(h.get_value())).points) + [h.get_value() * UP, ORIGIN]
            trsframe.set_points_as_corners(tarcpts)
            trsframe.rotate_about_origin(PI / 2, axis = RIGHT)
            tlsframe = trsframe.copy().rotate_about_origin(ang.get_value(), axis = OUT)
            temphemframe = VGroup(
                tlsframe, trsframe,
                ParametricFunction(lambda t: (1 - h.get_value() ** 2) ** 0.5 * np.array([np.cos(t), np.sin(t), 0]) + h.get_value() * OUT, t_min = 0, t_max = ang.get_value()),
            ).set_color(BLUE)
            tempthem.shift(h.get_value() * OUT)
            tvhem = VGroup(temphem, tempbhem, tempthem, temphemframe)
            tvhem.scale(2)
            tvhem.shift(3 * RIGHT)
            obj.become(tvhem)
        vhem.add_updater(hemupd)
        self.play(Write(axes))
        self.play(Write(punccone), Write(vhem))
        self.move_camera(phi = PI / 3)
        self.play(vhem.shift, 3 * RIGHT, punccone.shift, 3 * LEFT)
        self.wait()
        self.add(vhem)
        self.add(punccone)
        self.play(ang.increment_value, -PI / 2, run_time = 5)
        #self.begin_ambient_camera_rotation(rate = 1)
        vhem.clear_updaters()
        punccone.clear_updaters()
        self.wait(5)

class Revolutionv(SpecialThreeDScene):
    def construct(self):
        dsk_chk_colors = [GREEN_E]
        axes = self.get_axes()
        h = ValueTracker(1)
        ang = 2 * PI
        cne = Cone(v_max = ang, u_max = h.get_value(), checkerboard_colors = [RED_D, RED_E])
        cne.set_stroke(width = 0)
        cyl = Cylinder(v_max = ang, u_max = h.get_value(), checkerboard_colors = [BLUE_E])
        cyl.set_stroke(width = 0)
        bcone = ParametricDisc(v_max = ang)
        tcone = ParametricDisc(u_min = h.get_value(), v_max = ang, checkerboard_colors = dsk_chk_colors)
        tcone.shift(h.get_value() * OUT)
        coneframe = VGroup(
            ParametricFunction(lambda t: np.array([np.cos(t), np.sin(t), h.get_value()]), t_min = 0, t_max = ang),
            ParametricFunction(lambda t: np.array([h.get_value() * np.cos(t), h.get_value() * np.sin(t), h.get_value()]), t_min = 0, t_max = ang),
        ).set_color(BLUE)
        punccone = VGroup(cne, cyl, bcone, tcone, coneframe)
        punccone.scale(2, about_point = ORIGIN)
        def ghost(obj):
            res = obj.copy()
            res.set_fill(BLUE_E, opacity = 0.125)
            res.set_stroke(WHITE, width = 0.5, opacity = 0.125)
            return res
        #cone_ghost = ghost(punccone)
        #cone_ghost.shift(3 * LEFT)
        def cneupd(obj):
            tcne, tcyl, tbcone, ttcone, tconeframe = obj
            tcne.become(Cone(v_max = ang, u_max = h.get_value(), checkerboard_colors = [RED_E]).set_stroke(width = 0))
            tcyl.become(Cylinder(v_max = ang, u_max = h.get_value()).set_stroke(width = 0))
            tbcone.become(ParametricDisc(v_max = ang))
            ttcone.become(ParametricDisc(u_min = h.get_value(), v_max = ang, checkerboard_colors = dsk_chk_colors))
            ttcone.shift(h.get_value() * OUT)
            tconeframe.become(
                VGroup(
                    ParametricFunction(lambda t: np.array([np.cos(t), np.sin(t), h.get_value()]), t_min = 0, t_max = ang),
                    ParametricFunction(lambda t: np.array([h.get_value() * np.cos(t), h.get_value() * np.sin(t), h.get_value()]), t_min = 0, t_max = ang),
                ).set_color(BLUE)
            )
            tpunccone = VGroup(tcne, tcyl, tbcone, ttcone, tconeframe)
            tpunccone.scale(2, about_point = ORIGIN)
            tpunccone.shift(3 * LEFT)
            obj.become(tpunccone)
        hem = HemiSphere(u_max = np.arcsin(h.get_value()), v_max = ang)
        hem.set_stroke(width = 0)
        bhem = ParametricDisc(v_max = ang)
        them = ParametricDisc(u_max = (1 - h.get_value() ** 2) ** 0.5, v_max = ang, checkerboard_colors = dsk_chk_colors)
        hemframe = VGroup(
            ParametricFunction(lambda t: (1 - h.get_value() ** 2) ** 0.5 * np.array([np.cos(t), np.sin(t), 0]) + h.get_value() * OUT, t_min = 0, t_max = ang),
        ).set_color(BLUE)
        them.shift(h.get_value() * OUT)
        vhem = VGroup(hem, bhem, them, hemframe)
        vhem.scale(2, about_point = ORIGIN)
        #hem_ghost = ghost(vhem)
        #hem_ghost.shift(3 * RIGHT)
        def hemupd(obj):
            temphem, tempbhem, tempthem, temphemframe = obj
            temphem.become(HemiSphere(u_max = np.arcsin(h.get_value()), v_max = ang).set_stroke(width = 0))
            tempbhem.become(ParametricDisc(v_max = ang))
            tempthem.become(ParametricDisc(u_max = (1 - h.get_value() ** 2) ** 0.5, v_max = ang, checkerboard_colors = dsk_chk_colors))
            temphemframe = VGroup(
                ParametricFunction(lambda t: (1 - h.get_value() ** 2) ** 0.5 * np.array([np.cos(t), np.sin(t), 0]) + h.get_value() * OUT, t_min = 0, t_max = ang),
            ).set_color(BLUE)
            tempthem.shift(h.get_value() * OUT)
            tvhem = VGroup(temphem, tempbhem, tempthem, temphemframe)
            tvhem.scale(2, about_point = ORIGIN)
            tvhem.shift(3 * RIGHT)
            obj.become(tvhem)
        self.play(Write(axes))
        self.begin_ambient_camera_rotation()
        self.play(Write(punccone), Write(vhem))
        self.wait()
        self.move_camera(phi = PI / 4)
        self.wait()
        self.play(vhem.shift, 3 * RIGHT, punccone.shift, 3 * LEFT)
        self.wait()
        #self.add(hem_ghost, cone_ghost)
        vhem.add_updater(hemupd)
        punccone.add_updater(cneupd)
        self.add(vhem)
        self.add(punccone)
        self.play(h.increment_value, -0.5, run_time = 5)
        self.stop_ambient_camera_rotation()
        vhem.clear_updaters()
        punccone.clear_updaters()
        self.wait(10)

class TestRevolution(Scene):
    def construct(self):
        h = 0.5
        a = Arc(angle = np.arcsin(h))
        test = VMobject(fill_color = BLUE, fill_opacity = 0.875)
        pts = [ORIGIN] + list(a.points) + [h * UP, ORIGIN]
        test.set_points_as_corners(pts)
        '''test.add(
            Line(ORIGIN, RIGHT),
            a,
            Line(a.get_end(), h * UP),
            Line(h * UP, ORIGIN)
        )'''
        test.scale(2)
        self.play(Write(test), run_time = 2)
        self.wait(5)

class Revolutionv2(SpecialThreeDScene):
    def construct(self):
        dsk_chk_colors = [GREEN_E]
        axes = self.get_axes()
        h = ValueTracker(1)
        ang = 2 * PI
        cyl = Cylinder(v_max = ang, u_max = 1, checkerboard_colors = [BLUE_D, BLUE_E])
        cne = Cone(v_max = ang, u_max = 1, checkerboard_colors = [RED_D, RED_E])
        cnebase = ParametricDisc()
        puncyl = VGroup(cyl, cne, cnebase)
        puncyl.set_fill(BLUE_E, opacity = 0.5)
        puncyl.set_stroke(width = 0.5, opacity = 0.5)
        puncyl.scale(2, about_point = ORIGIN)
        puncyl.shift(3 * LEFT)
        hemi = HemiSphere(u_max = PI / 2, v_max = ang, checkerboard_colors = [BLUE_D, BLUE_E])
        hembase = ParametricDisc()
        hem = VGroup(hemi, hembase)
        hem.set_fill(BLUE_E, opacity = 0.5)
        hem.set_stroke(width = 0.5, opacity = 0.5)
        hem.scale(2, about_point = ORIGIN)
        hem.shift(3 * RIGHT)
        def conef(j):
            return VGroup(
            ParametricDisc(u_min = j, v_max = ang, checkerboard_colors = dsk_chk_colors).set_fill(GREEN_E).shift(j * OUT),
            ParametricFunction(lambda t: np.array([np.cos(t), np.sin(t), j]), t_min = 0, t_max = ang).set_color(BLUE),
            ParametricFunction(lambda t: np.array([j * np.cos(t), j * np.sin(t), j]), t_min = 0, t_max = ang).set_color(BLUE),).scale(2, about_point = ORIGIN).shift(3 * LEFT)
        coneframe = conef(h.get_value())
        def coneupd(obj):
            obj.become(conef(h.get_value()))
        coneframe.add_updater(coneupd)
        def hemf(j):
            return VGroup(
                ParametricDisc(u_max = (1 - j ** 2) ** 0.5, v_max = ang, checkerboard_colors = dsk_chk_colors).shift(j * OUT),
                ParametricFunction(lambda t: (1 - j ** 2) ** 0.5 * np.array([np.cos(t), np.sin(t), 0]) + j * OUT, t_min = 0, t_max = ang).set_color(BLUE),
            ).scale(2, about_point = ORIGIN).shift(3 * RIGHT)
        hemframe = hemf(h.get_value())
        def hemupd(obj):
            obj.become(hemf(h.get_value()))
        hemframe.add_updater(hemupd)
        self.move_camera(phi = 70 * DEGREES)
        self.play(Write(puncyl), Write(hem))
        self.add(coneframe, hemframe)
        self.wait()
        self.play(h.increment_value, -0.5, run_time = 5)
        self.wait(5)

class Revolutionv3(SpecialThreeDScene):
    def construct(self):
        oopac = 0.125
        aopac = 0.25
        axes = self.get_axes()
        h = ValueTracker(1)
        ang = 2 * PI
        Cyl = Cylinder(v_max = ang, u_max = 1)
        Cne = Cone(v_max = ang, u_max = 1)
        BseCne = ParametricDisc()
        PuncCone = VGroup(Cyl, Cne, BseCne)
        PuncCone.set_fill(BLUE_E, opacity = oopac)
        PuncCone.set_stroke(WHITE, width = 0.25, opacity = oopac)
        PuncCone.scale(2, about_point = ORIGIN)
        PuncCone.shift(3 * LEFT)
        Hem = HemiSphere(u_max = PI / 2, v_max = ang)
        HemBse = ParametricDisc()
        HemGrp = VGroup(Hem, HemBse)
        HemGrp.set_fill(BLUE_E, opacity = oopac)
        HemGrp.set_stroke(WHITE, width = 0.25, opacity = oopac)
        HemGrp.scale(2, about_point = ORIGIN)
        HemGrp.shift(3 * RIGHT)
        def makecone(j):
            cne = Cone(v_max = ang, u_max = j)
            cne.set_fill(BLUE_E, opacity = aopac)
            cyl = Cylinder(v_max = ang, u_max = j)
            cyl.set_fill(BLUE_E, opacity = aopac)
            bcone = ParametricDisc()
            bcone.set_fill(BLUE_E, opacity = aopac)
            tcone = ParametricDisc(u_min = j, v_max = ang)
            tcone.shift(j * OUT)
            tcone.set_fill(GREEN, opacity = aopac)
            res = VGroup(cne, cyl, bcone, tcone)
            res.set_stroke(width = 0)
            res.scale(2, about_point = ORIGIN)
            res.shift(3 * LEFT)
            return res
        MadeCone = makecone(h.get_value())
        def makesphere(j):
            hem = HemiSphere(u_max = np.arcsin(j), v_max = ang)
            hem.set_fill(BLUE_E, opacity = aopac)
            bhem = ParametricDisc()
            bhem.set_fill(BLUE_E, opacity = aopac)
            them = ParametricDisc(u_max = (1 - j ** 2) ** 0.5, v_max = ang)
            them.shift(j * OUT)
            them.set_fill(GREEN, opacity = aopac)
            res = VGroup(hem, bhem, them)
            res.set_stroke(width = 0)
            res.scale(2, about_point = ORIGIN)
            res.shift(3 * RIGHT)
            return res
        MadeSphere = makesphere(h.get_value())
        self.move_camera(phi = 75 * DEGREES)
        self.play(Write(PuncCone), Write(HemGrp))
        def coneupd(obj):
            obj.become(makecone(h.get_value()))
        def sphupd(obj):
            obj.become(makesphere(h.get_value()))
        MadeCone.add_updater(coneupd)
        MadeSphere.add_updater(sphupd)
        self.add(MadeCone, MadeSphere)
        #self.play(PuncCone.set_fill, {"opacity": 0.125}, HemGrp.set_fill, {"opacity": 0.125})
        self.play(h.increment_value, -0.75, run_time = 10)
        MadeCone.clear_updaters()
        MadeSphere.clear_updaters()
        self.wait(5)

class Revolutionv4(SpecialThreeDScene):
    def construct(self):
        oopac = 0.125
        aopac = 0.375
        #axes = self.get_axes()
        h = ValueTracker(1)
        ang = 2 * PI
        Cyl = Cylinder(v_max = ang, u_max = 1)
        Cne = Cone(v_max = ang, u_max = 1)
        BseCne = ParametricDisc()
        PuncCone = VGroup(Cyl, Cne, BseCne)
        PuncCone.set_fill(BLUE_E, opacity = oopac)
        PuncCone.set_stroke(WHITE, width = 0.25, opacity = oopac)
        PuncCone.scale(2, about_point = ORIGIN)
        PuncCone.shift(3 * LEFT)
        Hem = HemiSphere(u_max = PI / 2, v_max = ang)
        HemBse = ParametricDisc()
        HemGrp = VGroup(Hem, HemBse)
        HemGrp.set_fill(BLUE_E, opacity = oopac)
        HemGrp.set_stroke(WHITE, width = 0.25, opacity = oopac)
        HemGrp.scale(2, about_point = ORIGIN)
        HemGrp.shift(3 * RIGHT)
        def makecone(j):
            cne = Cone(v_max = ang, u_max = j)
            cne.set_fill(BLUE_E, opacity = aopac)
            cyl = Cylinder(v_max = ang, u_max = j)
            cyl.set_fill(BLUE_E, opacity = aopac)
            bcone = ParametricDisc()
            bcone.set_fill(BLUE_E, opacity = aopac)
            tcone = ParametricDisc(u_min = j, v_max = ang)
            tcone.shift(j * OUT)
            tcone.set_fill(GREEN, opacity = aopac)
            res = VGroup(cne, cyl, bcone, tcone)
            res.set_stroke(width = 0)
            res.scale(2, about_point = ORIGIN)
            res.shift(3 * LEFT)
            return res
        MadeCone = makecone(h.get_value())
        def makesphere(j):
            hem = HemiSphere(u_max = np.arcsin(j), v_max = ang)
            hem.set_fill(BLUE_E, opacity = aopac)
            bhem = ParametricDisc()
            bhem.set_fill(BLUE_E, opacity = aopac)
            them = ParametricDisc(u_max = (1 - j ** 2) ** 0.5, v_max = ang)
            them.shift(j * OUT)
            them.set_fill(GREEN, opacity = aopac)
            res = VGroup(hem, bhem, them)
            res.set_stroke(width = 0)
            res.scale(2, about_point = ORIGIN)
            res.shift(3 * RIGHT)
            return res
        MadeSphere = makesphere(h.get_value())
        self.move_camera(phi = 75 * DEGREES)
        #self.play(Write(PuncCone), Write(HemGrp))
        def coneupd(obj):
            obj.become(makecone(h.get_value()))
        def sphupd(obj):
            obj.become(makesphere(h.get_value()))
        self.play(Write(MadeCone), Write(MadeSphere), Write(PuncCone), Write(HemGrp))
        MadeCone.add_updater(coneupd)
        MadeSphere.add_updater(sphupd)
        self.add(MadeCone, MadeSphere)
        #self.play(PuncCone.set_fill, {"opacity": 0.125}, HemGrp.set_fill, {"opacity": 0.125})
        self.play(h.increment_value, -0.75, run_time = 10)
        MadeCone.clear_updaters()
        MadeSphere.clear_updaters()
        self.wait(5)

class CavExample(SpecialThreeDScene):
    def construct(self):
        #axes = ThreeDAxes()
        def corprism(j):
            res = VMobject(fill_opacity = 0.75, fill_color = BLUE, storke_width = 0.0)
            j = 1 if j > 1 else j
            t = 1 - j
            face1 = res.copy()
            face1.set_points_as_corners([ORIGIN, RIGHT, RIGHT + UP, UP, ORIGIN])
            face2 = res.copy()
            face2.set_points_as_corners([ORIGIN, t * RIGHT, t * RIGHT + t * UP, t * UP, ORIGIN])
            face2.shift(j * OUT)
            face2.set_fill(GREEN)
            face3 = res.copy()
            face3.set_points_as_corners([ORIGIN, RIGHT, t * RIGHT + j * OUT, j * OUT, ORIGIN])
            face4 = res.copy()
            face4.set_points_as_corners([RIGHT, t * RIGHT + j * OUT, t * RIGHT + t * UP + j * OUT, RIGHT + UP, RIGHT])
            face5 = res.copy()
            face5.set_points_as_corners([RIGHT + UP, t * RIGHT + t * UP + j * OUT, t * UP + j * OUT, UP, RIGHT + UP])
            face6 = res.copy()
            face6.set_points_as_corners([ORIGIN, j * OUT, j * OUT + t * UP, UP, ORIGIN])
            resobj = VGroup(face1, face2, face3, face4, face5, face6)
            resobj.set_stroke(WHITE, width = 0.5)
            return resobj
        def cenprism(j):
            res = VMobject(fill_opacity = 0.75, fill_color = BLUE, storke_width = 0.0)
            j = 1 if j > 1 else j
            t = 1 - j
            face1 = res.copy()
            face1.set_points_as_corners([ORIGIN, RIGHT, RIGHT + UP, UP, ORIGIN])
            face2 = res.copy()
            face2.set_points_as_corners([ORIGIN, t * RIGHT, t * RIGHT + t * UP, t * UP, ORIGIN])
            face2.shift((1 - t) * RIGHT / 2 + (1 - t) * UP / 2 + j * OUT)
            face2.set_fill(GREEN)
            face3 = res.copy()
            face3.set_points_as_corners([ORIGIN, RIGHT, RIGHT * (1 + t) / 2 + UP * (1 - t) / 2 + j * OUT, RIGHT * (1 - t) / 2 + UP * (1 - t) / 2 + j * OUT, ORIGIN])
            face4 = res.copy()
            face4.set_points_as_corners([RIGHT, RIGHT * (1 + t) / 2 + UP * (1 - t) / 2 + j * OUT, RIGHT * (1 + t) / 2 + UP * (1 + t) / 2 + j * OUT, RIGHT + UP, RIGHT])
            face5 = res.copy()
            face5.set_points_as_corners([RIGHT + UP, RIGHT * (1 + t) / 2 + UP * (1 + t) / 2 + j * OUT, RIGHT * (1 - t) / 2 + UP * (1 + t) / 2 + j * OUT, UP, RIGHT + UP])
            face6 = res.copy()
            face6.set_points_as_corners([UP, RIGHT * (1 - t) / 2 + UP * (1 + t) / 2 + j * OUT, RIGHT * (1 - t) / 2 + UP * (1 - t) / 2 + j * OUT, ORIGIN, UP])
            resobj = VGroup(face1, face2, face3, face4, face5, face6)
            resobj.set_stroke(WHITE, width = 0.5)
            return resobj
        h = ValueTracker(1)
        e = cenprism(h.get_value())
        e.scale(3, about_point = ORIGIN)
        e.shift(4.5 * LEFT + 1.5 * DOWN)
        r = corprism(h.get_value())
        r.scale(3, about_point = ORIGIN)
        r.shift(1.5 * RIGHT + 1.5 * DOWN)
        eg = cenprism(1).set_fill(opacity = 0.0625).scale(3, about_point = ORIGIN).shift(4.5 * LEFT + 1.5 * DOWN)
        rg = corprism(1).set_fill(opacity = 0.0625).scale(3, about_point = ORIGIN).shift(1.5 * RIGHT + 1.5 * DOWN)
        def cenupd(obj):
            res = cenprism(h.get_value())
            res.scale(3, about_point = ORIGIN).shift(4.5 * LEFT + 1.5 * DOWN)
            obj.become(res)
        def corupd(obj):
            res = corprism(h.get_value())
            res.scale(3, about_point = ORIGIN).shift(1.5 * RIGHT + 1.5 * DOWN)
            obj.become(res)
        e.add_updater(cenupd)
        r.add_updater(corupd)
        self.move_camera(phi = 75 * DEGREES)
        self.begin_ambient_camera_rotation(rate = 0.125)
        self.play(
            #Write(axes),
            Write(rg), Write(eg), Write(e), Write(r)
        )
        self.add(e, r)
        self.wait()
        self.play(h.increment_value, -0.75, run_time = 10)
        self.play(h.increment_value, 0.75, run_time = 10)
        self.play(h.increment_value, -0.75, run_time = 10)
        self.wait(30)

class Revolutionv5(SpecialThreeDScene):
    def construct(self):
        axes = self.get_axes()
        self.add(axes)
        fullcurve = ParametricFunction(lambda t: np.array([t, 1 / t, 0]), t_min = 0.01, t_max = 100).set_color(BLUE)
        horncurve = VGroup(
            ParametricFunction(lambda t: np.array([t, 1, 0]), t_min = 0, t_max = 1),
            ParametricFunction(lambda t: np.array([t, 1 / t, 0]), t_min = 1, t_max = 100),
        ).set_color(YELLOW)
        self.play(Write(fullcurve))
        self.wait()
        self.play(FadeOut(fullcurve), FadeIn(horncurve))
        self.wait()
        h = GabHorn()
        c = Cylinder(checkerboard_colors = [BLUE, GREEN, YELLOW, ORANGE], resolution = (12, 24))
        c.set_stroke(width = 0, opacity = 0.25)
        c.rotate_about_origin(PI / 2, axis = UP)
        horn = VGroup(h, c)
        horn.set_fill([BLUE_E, GREEN_E], opacity = 0.5)
        self.play(Write(horn))
        self.wait()
        self.begin_ambient_camera_rotation()
        self.move_camera(
            **self.default_angled_camera_position,
            run_time = 2,
        )
        self.play(
            horn.rotate_about_origin, PI / 2, {"axis": RIGHT},
            horncurve.rotate_about_origin, PI / 2, {"axis": RIGHT},
        )
        self.wait()
        cy = Cylinder(u_max = 1 / np.sqrt(2))
        tbase = ParametricDisc().shift(OUT / np.sqrt(2))
        bbase = ParametricDisc()
        tbase.set_shade_in_3d(True)
        bbase.set_shade_in_3d(True)
        cyl = VGroup(cy, tbase, bbase)
        cyl.scale(np.sqrt(2), about_point = ORIGIN)
        #cyl.rotate_about_origin(-PI / 2, axis = RIGHT)
        cyl.shift(3 * LEFT)
        cyl.set_stroke(width = 0.25, opacity = 0.5)
        cyl.set_fill([BLUE_D, BLUE_E], opacity = 0.25)
        #self.move_camera(phi = -45 * DEGREES, theta = -90 * DEGREES)
        #self.set_to_default_angled_camera_orientation()
        self.play(Write(cyl))
        self.wait()
        self.stop_ambient_camera_rotation()
        self.move_camera(
            phi = 90 * DEGREES,
            theta = -90 * DEGREES,
        )
        tworoottwo = TexMobject("2\\sqrt{2}").rotate(90 * DEGREES, RIGHT)
        onetext = TexMobject("1").rotate(90 * DEGREES, RIGHT)
        twolabelline = Line(OUT + (3 - np.sqrt(2)) * LEFT, OUT + (3 + np.sqrt(2)) * LEFT).set_color(YELLOW)
        onelabelline = Line((3 + np.sqrt(2)) * LEFT, (3 + np.sqrt(2)) * LEFT + OUT).set_color(RED)
        onetext.next_to(onelabelline, LEFT)
        tworoottwo.next_to(twolabelline, OUT)
        self.play(
            *[FadeIn(obj) for obj in [onetext, tworoottwo, onelabelline, twolabelline]]
        )
        self.wait()
        def cavgrp(j):
            x = 1 / j
            tline = Line(j * UP, x * RIGHT + j * UP).set_color(GOLD_E)
            tline.rotate_about_origin(PI / 2, axis = RIGHT)
            tcyl = Cylinder(radius = j, u_max = x)
            tcyl.set_stroke(width = 0.125, opacity = 0.25)
            tcyl.set_fill(GOLD_D, opacity = 0.25)
            tcyl.rotate_about_origin(PI / 2, axis = UP)
            tdisc = ParametricDisc().scale(np.sqrt(2), about_point = ORIGIN).shift(j * OUT)
            #tdisc.rotate_about_origin(-PI / 2, axis = RIGHT)
            tdisc.shift(3 * LEFT)
            tplane = FlatPlane(u_min = 0, v_min = -PI * j, u_max = x, v_max = PI * j, checkerboard_colors = [GOLD_C, GOLD_B, GOLD_A]).shift(j * OUT)
            #tplane.rotate_about_origin(-PI / 2, axis = RIGHT)
            tplane.set_fill(opacity = 0.5)
            res = VGroup(tline, tcyl, tdisc, tplane)
            return res
        def grpupd(obj):
            obj.become(cavgrp(y.get_value()))
        y = ValueTracker(0.75)
        #self.move_camera(phi = 10 * DEGREES, theta = 10 * DEGREES, frame_center = ORIGIN)
        #self.play(ApplyMethod(self.camera.phi_tracker.increment_value, 10 * DEGREES))
        '''self.play(
            self.camera.theta_tracker.increment_value, -0 * DEGREES,
            self.camera.phi_tracker.increment_value, -45 * DEGREES,
            )'''
        #self.set_to_default_angled_camera_orientation()
        mgp = cavgrp(y.get_value())
        tempcyl = mgp[1].copy()
        tempplane = mgp[3].copy()
        self.play(Write(mgp[0]))
        self.wait()
        coordlabel = TexMobject("(h^{-1},h)").scale(0.75).rotate(90 * DEGREES, RIGHT).next_to(mgp[0], RIGHT, buff = 0.125)
        self.play(Write(coordlabel))
        self.wait()
        self.begin_ambient_camera_rotation()
        self.move_camera(
            **self.default_angled_camera_position,
            run_time = 2,
        )
        lhorn = VGroup(
            GabHorn(v_max = PI).rotate_about_origin(-PI / 2, axis = RIGHT),
            Cylinder(checkerboard_colors = [BLUE, GREEN, YELLOW, ORANGE], v_max = PI, resolution = (12, 24)).set_stroke(width = 0, opacity = 0.25).rotate_about_origin(PI / 2, axis = UP)
        ).set_fill([BLUE_E, GREEN_E], opacity = 0.5)
        rhorn = VGroup(
            GabHorn(v_min = PI, v_max = 2 * PI).rotate_about_origin(-PI / 2, axis = RIGHT),
            Cylinder(checkerboard_colors = [BLUE, GREEN, YELLOW, ORANGE], v_min = PI, v_max = 2 * PI, resolution = (12, 24)).set_stroke(width = 0, opacity = 0.25).rotate_about_origin(PI / 2, axis = UP)
        ).set_fill([BLUE_E, GREEN_E], opacity = 0.5)
        self.add(lhorn, rhorn)
        self.remove(horn)
        #self.play(rhorn.rotate, PI / 2, {"axis": RIGHT, "about_point": OUT})
        self.play(
            Rotate(
                rhorn, -PI,
                axis = RIGHT,
                about_point = OUT,
                run_time = 2
            ),
            VFadeOut(rhorn, run_time = 2)
        )
        self.wait()
        cylang = 0
        yval = y.get_value()
        startang = PI
        angval = ValueTracker(PI)
        def cylcreator(j):
            obj = Cylinder(radius = yval, u_max = 1 / yval, v_min = startang, v_max = j)
            obj.set_stroke(width = 0.125, opacity = 0.25)
            obj.set_fill(GOLD_D, opacity = 0.25)
            obj.rotate_about_origin(PI / 2, axis = UP)
            return obj
        def cylupd(obj):
            obj.become(cylcreator(angval.get_value()))
        rotcyl = cylcreator(angval.get_value())
        rotcyl.add_updater(cylupd)
        self.add(rotcyl)
        self.play(angval.increment_value, 2 * PI, run_time = 2)
        self.wait()
        self.add(tempcyl)
        self.remove(rotcyl)
        rotcyl.clear_updaters()
        self.play(ReplacementTransform(tempcyl.copy(), tempplane))
        self.wait()
        self.stop_ambient_camera_rotation()
        self.move_camera(
            phi = 0 * DEGREES,
            theta = -90 * DEGREES,
        )
        lenlabel = TexMobject("1/h")
        #brelabel = TexMobject("\\frac{2\\pi}{h}")
        brelabel = TexMobject("2\\pi h").rotate_about_origin(PI / 2)
        lenbrace = Brace(Line(PI * UP * yval + yval * OUT, PI * UP * yval + RIGHT / yval + yval * OUT), UP, buff = 0.125)
        brebrace = Brace(Line(PI * UP * yval + yval * OUT, PI * DOWN * yval + yval * OUT), LEFT, buff = 0.125)
        lenlabel.next_to(lenbrace, UP, buff = 0.125)
        brelabel.next_to(brebrace, LEFT, buff = 0.125)
        self.play(*[Write(obj) for obj in [lenlabel, brelabel, lenbrace, brebrace]])
        self.wait()
        self.begin_ambient_camera_rotation()
        self.move_camera(
            **self.default_angled_camera_position,
            run_time = 2,
        )
        self.add(mgp[1], mgp[3])
        self.remove(tempcyl, tempplane)
        self.play(Write(mgp[2]))
        self.wait()
        self.play(*[FadeOut(obj) for obj in [coordlabel, onetext, tworoottwo, onelabelline, twolabelline, lenlabel, brelabel, lenbrace, brebrace]])
        self.wait()
        mgp.add_updater(grpupd)
        self.add(mgp)
        self.play(y.increment_value, -0.65, run_time = 5, rate_func = linear)
        self.wait()
        self.play(y.increment_value, 0.65, run_time = 5, rate_func = linear)
        self.stop_ambient_camera_rotation()
        self.move_camera(
            phi = 90 * DEGREES,
            theta = -90 * DEGREES,
            run_time = 2
        )
        mgp.clear_updaters()
        volumetext = TextMobject("Volume$=2\\pi$").rotate(PI / 2, axis = RIGHT).move_to(3 * RIGHT + 2 * OUT)
        self.play(Write(volumetext))
        self.wait(5)

class Revolutionv6(SpecialThreeDScene):
    def construct(self):
        axes = self.get_axes()
        self.add(axes)
        torri = ImageMobject("Torri").scale(2.5)
        torriname = TextMobject("Evangelista Torricelli")
        torri.shift(3 * LEFT)
        torriname.next_to(torri, DOWN)
        fullcurve = ParametricFunction(lambda t: np.array([t, 1 / t, 0]), t_min = 0.01, t_max = 100).set_color(BLUE)
        self.play(Write(fullcurve))
        self.wait()
        curveeqn = TexMobject("y=\\frac{1}{x}").move_to(3 * UR)
        self.play(Write(curveeqn))
        self.wait()
        h = GabHorn()
        self.play(Write(h))
        self.wait()
        self.play(Write(torriname), FadeIn(torri))
        self.add_fixed_in_frame_mobjects(torri, torriname)
        self.begin_ambient_camera_rotation(rate = 0.02)
        self.move_camera(
            **self.default_angled_camera_position,
            run_time = 2,
        )
        self.play(h.set_fill, {"opacity": 0.75})
        self.wait(30)
        self.stop_ambient_camera_rotation()
        self.wait()

if __name__ == "__main__":
    module_name = os.path.basename(__file__)
    #command_A = "manim -p  --video_dir ~/Downloads/  "
    clear_cmd = "cls"
    command_A = "python -m manim "
    command_B = module_name + " " + "CavExample" + " -p"
    #command_B = module_name + " " + "IntroductionScene" + " -p"
    #command_B = module_name + " -p"
    os.system(clear_cmd)
    os.system(command_A + command_B)