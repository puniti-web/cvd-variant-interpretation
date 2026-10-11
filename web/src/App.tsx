import { BrowserRouter, Route, Routes } from 'react-router-dom'
import { Shell } from './components/Shell.tsx'
import { AboutPage } from './features/about/AboutPage.tsx'
import { BenchPage } from './features/bench/BenchPage.tsx'
import { LandingPage } from './features/landing/LandingPage.tsx'
import { ExplorePage } from './features/explore/ExplorePage.tsx'
import { ModelsPage } from './features/models/ModelsPage.tsx'
import { ReportPage } from './features/report/ReportPage.tsx'
import { StructurePage } from './features/structure/StructurePage.tsx'
import { EmptyState } from './components/ui.tsx'

export default function App() {
  return (
    <BrowserRouter>
      <Routes>
        <Route path="/" element={<LandingPage />} />
        <Route element={<Shell />}>
          <Route path="bench" element={<BenchPage />} />
          <Route path="explore" element={<ExplorePage />} />
          <Route path="structure/:id" element={<StructurePage />} />
          <Route path="report/:id" element={<ReportPage />} />
          <Route path="models" element={<ModelsPage />} />
          <Route path="about" element={<AboutPage />} />
          <Route path="*" element={<EmptyState title="This drawer is empty" detail="The address does not match a page." />} />
        </Route>
      </Routes>
    </BrowserRouter>
  )
}
